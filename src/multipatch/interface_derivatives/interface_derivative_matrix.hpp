// SPDX-License-Identifier: MIT

#pragma once

#include <Eigen>

#include <ddc/ddc.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "matching_idx_slice.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
// #include "single_interface_derivatives_calculator.hpp"
#include "types.hpp"

/**
 * TODO: 
 *  - How does it work? 
 *      - Constructor
 *          - Take a list of "parallel" interfaces. 
 *          - For each "conforming" lines, get the coefficients a and b.
 *          - Fill in the matrix. 
 *      - operator()
 *          - get the coefficient c. 
 *          - Fill in the vector. 
 *          - Inverse the matrix system. (See how it is done in the SplineBuilder.)
 * 
 *  - Remarks: 
 *      - Should we restrain to the conforming case first? 
 *      - And then the non-conforming case in another operator? 
 *      - Treat only "one block line/direction"? 
 *      - And build another operator to manage the whole geometry? 
 *      - Including treatment for T-joints? 
 */

template <
        class Connectivity,
        class Grid1D,
        template <typename P>
        typename ValuesOnPatch,
        template <typename P>
        typename DerivsOnPatch,
        ddc::BoundCond LowerBound = ddc::BoundCond::HERMITE,
        ddc::BoundCond UpperBound = ddc::BoundCond::HERMITE,
        class ExecSpace = Kokkos::DefaultHostExecutionSpace,
        class... Patches>
class InterfaceDerivativeMatrix
{
    // using interface_collection = typename Connectivity::interface_collection; // TypeSeq
    using interface_collection =
            typename Connectivity::get_all_interfaces_along_direction_t<Grid1D>;
    using all_patches = typename Connectivity::all_patches; // TypeSeq
    using grid_collection = collect_grids_on_dim_t<Grid1D>;

    // TODO: remove Interfaces with OutsideEdge
    using inner_interface_collection = ddc::detail::TypeSeq<
            ddc::type_seq_element_t<8, interface_collection>,
            ddc::type_seq_element_t<9, interface_collection>>;


    // HELPFUL ALIASES ===========================================================================
    template <typename Grid>
    using get_patch_on_grid = find_patch_t<Grid, all_patches>;


    template <typename T>
    struct get_grid;

    template <template <typename P> typename T, typename Grid>
    struct get_grid<T<Grid>>
    {
        using type = Grid;
    };

    template <typename T>
    using get_grid_t = typename get_grid<T>::type;


    template <typename IdxRange1D>
    using get_patch_of_idx_range_1d = get_patch_on_grid<get_grid_t<IdxRange1D>>;


    template <template <typename P> typename T, class InterfacesTypeSeq>
    struct get_tuple_on_interfaces;

    template <template <typename P> typename T, class InterfacesTypeSeq>
    using get_tuple_on_interfaces_t = typename get_tuple_on_interfaces<T, InterfacesTypeSeq>::type;

    template <template <typename P> typename T, class... Interfaces>
    struct get_tuple_on_interfaces<T, ddc::detail::TypeSeq<Interfaces...>>
    {
        using type = std::tuple<T<Interfaces> const&...>;
    };

    template <typename InterfaceCollection>
    struct get_tuple_deriv_calculator;

    template <class... Interfaces>
    struct get_tuple_deriv_calculator<ddc::detail::TypeSeq<Interfaces...>>
    {
        using InterfaceCollection = ddc::detail::TypeSeq<Interfaces...>;
        using type = std::tuple<get_deriv_calculator_t<Interfaces, InterfaceCollection> const&...>;
    };

    template <class... Interfaces>
    using get_tuple_deriv_calculator_t
            = get_tuple_deriv_calculator<ddc::detail::TypeSeq<Interfaces...>>::type;

    template <class Interface, class InterfaceCollection>
    struct get_deriv_calculator
    {
        std::size_t nb_interface = ddc::type_seq_size_v<InterfaceCollection>;
        if constexpr (std::is_same_v<Interface, ddc::type_seq_element_t<0, InterfaceCollection>>) {
            using type = SingleInterfaceDerivativesCalculator<
                    Interface,
                    LowerBound,
                    ddc::BoundCond::HERMITE>;
        } else if (std::is_same_v<
                           Interface,
                           ddc::type_seq_element_t<nb_interface - 1, InterfaceCollection>>) {
            using type = SingleInterfaceDerivativesCalculator<
                    Interface,
                    ddc::BoundCond::HERMITE,
                    UpperBound>;
        } else {
            using type = SingleInterfaceDerivativesCalculator<Interface>;
        }
    }



    // ===========================================================================================

    // static constexpr std::size_t n_interfaces = ddc::type_seq_size_v<interface_collection>;
    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;
    static constexpr std::size_t n_values = n_inner_interfaces * 3 - 2;

    // using Matrix = gko::matrix::Csr<double>;
    using Matrix = gko::matrix::Dense<double>;



private:
    // Matrix m_matrix;
    std::shared_ptr<Matrix> m_matrix;
    std::shared_ptr<Matrix> m_vector;

    // std::array<std::array<double, n_inner_interfaces>, n_inner_interfaces> m_matrix;
    // std::array<double, n_inner_interfaces> m_vector;

    MultipatchType<IdxRangeOnPatch, Patches...> const& m_idx_ranges;



    // store it to not have to re-compute some coefficients.
    // using SingleInterfaceDerivativesCalculatorTuple = get_tuple_on_interfaces_t<
    //         SingleInterfaceDerivativesCalculator,
    //         inner_interface_collection>;
    // using SingleInterfaceDerivativesCalculatorTuple = std::tuple<
    //         SingleInterfaceDerivativesCalculator<
    //                 ddc::type_seq_element_t<0, inner_interface_collection>,
    //                 LowerBound,
    //                 ddc::BoundCond::HERMITE> const&,
    //         // SingleInterfaceDerivativesCalculator<
    //         //         ddc::type_seq_element_t<0, inner_interface_collection>> const&,
    //         SingleInterfaceDerivativesCalculator<
    //                 ddc::type_seq_element_t<n_inner_interfaces - 1, inner_interface_collection>,
    //                 ddc::BoundCond::HERMITE,
    //                 UpperBound> const&>;

    using SingleInterfaceDerivativesCalculatorTuple
            = get_tuple_deriv_calculator_t<inner_interface_collection>;
    SingleInterfaceDerivativesCalculatorTuple const& m_derivatives_calculators;

public:
    ~InterfaceDerivativeMatrix() {}

    InterfaceDerivativeMatrix(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            SingleInterfaceDerivativesCalculatorTuple const& derivatives_calculators)
        : m_idx_ranges(idx_ranges)
        , m_derivatives_calculators(derivatives_calculators)
    {
        set_matrix(std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
    }


    void solve(
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max)
    {
        std::tuple sorted_idx_ranges_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);
        using FirstIdxRange = std::tuple_element_t<0, decltype(sorted_idx_ranges_tuple)>;
        using FirstPatch = get_patch_of_idx_range_1d<FirstIdxRange>;
        // using FirstPatch = connectivity_details::FindPatch<Grid1D, all_patches>;
        auto par_idx_range = ddc::remove_dims_of(
                m_idx_ranges.template get<FirstPatch>(),
                std::get<0>(sorted_idx_ranges_tuple));

        ddc::for_each(par_idx_range, [&](auto const& idx) {
            std::array<double, n_inner_interfaces> derivs_at_interfaces
                    = this->solve_single_line(idx, function_values, derivs_min, derivs_max);

            update_derivatives(
                    derivs_min,
                    derivs_max,
                    derivs_at_interfaces,
                    function_values,
                    idx,
                    std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
        });
    }


private:
    template <std::size_t... I>
    void set_matrix(std::integer_sequence<std::size_t, I...>)
    {
        // Initialise to zero
        // for (int n(0); n < n_inner_interfaces; n++) {
        //     for (int m(0); m < n_inner_interfaces; m++) {
        //         m_matrix[n][m] = 0;
        //     }
        // }
        (set_line_matrix<I>(), ...);
    }

    // Fill in the line of the matrix corresponding to the given interface.
    template <std::size_t I>
    void set_line_matrix()
    {
        // auto [vals, cols, rows] = m_matrix->get_batch_csr();

        double coeff_right = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
        double coeff_left = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();
        const double coefs[] = {-coeff_left, 1, -coeff_right};

        // if (I == 0) {
        //     m_matrix[I][I + 1] = -coeff_right;
        //     m_matrix[I][I] = 1;
        // } else if (1 < I && I < n_inner_interfaces - 1) {
        //     m_matrix[I][I + 1] = -coeff_right;
        //     m_matrix[I][I] = 1;
        //     m_matrix[I][I - 1] = -coeff_left;
        // } else if (I == n_inner_interfaces - 1) {
        //     m_matrix[I][I - 1] = -coeff_left;
        //     m_matrix[I][I] = 1;
        // }
        for (auto dofs : {-1, 0, 1}) {
            if (0 <= I + dofs && I + dofs < n_inner_interfaces) {
                // vals(3 * I + dofs) = coefs[dofs + 1];
                // cols(3 * I + dofs) = I + dofs;
                // rows(I) = I - 1;
                matrix->at(I, I + dofs) = coefs[dofs + 1];
            }
        }
    }



    template <class OGrid1D>
    std::array<double, n_inner_interfaces> solve_single_line(
            Idx<OGrid1D> const& slice_idx,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max)
    {
        set_vector(
                slice_idx,
                function_values,
                derivs_min,
                derivs_max,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});

        // S = (I - M)^{-1} C = m_matrix^{-1} m_vector
        std::array<double, n_inner_interfaces> derivs_at_interfaces;
        // ...
        std::array<std::array<double, n_inner_interfaces>, n_inner_interfaces> inverse_matrix;
        double det = m_matrix[0][0] * m_matrix[1][1] - m_matrix[1][0] * m_matrix[0][1];
        inverse_matrix[0][0] = m_matrix[1][1] / det;
        inverse_matrix[0][1] = -m_matrix[0][1] / det;
        inverse_matrix[1][0] = -m_matrix[1][0] / det;
        inverse_matrix[1][1] = m_matrix[0][0] / det;

        derivs_at_interfaces[0]
                = inverse_matrix[0][0] * m_vector[0] + inverse_matrix[0][1] * m_vector[1];
        derivs_at_interfaces[1]
                = inverse_matrix[1][0] * m_vector[0] + inverse_matrix[1][1] * m_vector[1];

        return derivs_at_interfaces;
    }



    template <std::size_t... I, class OGrid1D>
    void set_vector(
            Idx<OGrid1D> const& slice_idx,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max,
            std::integer_sequence<std::size_t, I...>)
    {
        std::tuple sorted_idx_ranges_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);
        using FirstIdxRange = std::tuple_element_t<0, decltype(sorted_idx_ranges_tuple)>;
        using FirstPatch = get_patch_of_idx_range_1d<FirstIdxRange>;
        static_assert(
                (std::is_same_v<OGrid1D, typename FirstPatch::Grid1>)
                || (std::is_same_v<OGrid1D, typename FirstPatch::Grid2>)
                || (std::is_same_v<OGrid1D, ddc::Deriv<typename FirstPatch::Dim1>>)
                || (std::is_same_v<OGrid1D, ddc::Deriv<typename FirstPatch::Dim2>>));

        IdxRange<OGrid1D> idx_range_1d_first_patch(m_idx_ranges.template get<FirstPatch>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        // for (int n(0); n < n_inner_interfaces; n++) {
        //     m_vector->get_values()[n] = 0;
        // }
        (set_line_vector<I>(slice_idx_value, function_values, derivs_min, derivs_max), ...);
    }



    template <std::size_t I>
    void set_line_vector(
            int& slice_idx_1_value,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max)
    {
        // hyp: all_patches and sorted_idx_ranges_tuple have same organisation.
        std::tuple sorted_idx_ranges_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);
        static_assert(
                std::tuple_size_v<decltype(sorted_idx_ranges_tuple)> == n_inner_interfaces + 1);

        using Patch1 = get_patch_of_idx_range_1d<
                std::tuple_element_t<I, decltype(sorted_idx_ranges_tuple)>>;
        using Patch2 = get_patch_of_idx_range_1d<
                std::tuple_element_t<I + 1, decltype(sorted_idx_ranges_tuple)>>;

        // Define all the different index ranges
        auto const idx_range_1d_1 = std::get<I>(sorted_idx_ranges_tuple);
        auto const idx_range_1d_2 = std::get<I + 1>(sorted_idx_ranges_tuple);


        // using interfaces_typeseq =
        //         typename Connectivity::template find_connections_t<Patch1, Patch2>;
        // static_assert(ddc::type_seq_size_v<interfaces_typeseq> == 1);
        // using Interface = typename ddc::type_seq_element_t<0, interfaces_typeseq>;
        // MatchingIdxSlice<Interface> matching_idx(other_idx_range_1d_1, other_idx_range_1d_2);

        // Get the 2D index ranges of each patch where the function is interpolated.
        auto const idx_range_2d_function_1 = get_idx_range(function_values.template get<Patch1>());
        auto const idx_range_2d_function_2 = get_idx_range(function_values.template get<Patch2>());

        // Get the 1D index ranges parallel to the Interface.
        auto const other_idx_range_1d_function_1
                = ddc::remove_dims_of(idx_range_2d_function_1, idx_range_1d_1);
        auto const other_idx_range_1d_function_2
                = ddc::remove_dims_of(idx_range_2d_function_2, idx_range_1d_2);


        auto slice_indexes = get_slice_indexes<
                I,
                Patch1,
                Patch2>(slice_idx_1_value, sorted_idx_ranges_tuple, function_values);
        auto slice_idx_1 = std::get<0>(slice_indexes);
        auto slice_idx_2 = std::get<1>(slice_indexes);

        // Get a MultipatchField of the function on a perpendicular slice.
        auto function_slice_1
                = get_function_on_idx_range(idx_range_1d_1, slice_idx_1, function_values);
        auto function_slice_2
                = get_function_on_idx_range(idx_range_1d_2, slice_idx_2, function_values);

        double fct_coeff = std::get<I>(m_derivatives_calculators)
                                   .get_function_coefficients(function_slice_1, function_slice_2);

        // SUGESTION =============================================================================
        auto const idx_range_perp_1 = std::get<I>(sorted_idx_ranges_tuple);
        auto const idx_range_perp_2 = std::get<I + 1>(sorted_idx_ranges_tuple);

        ValuesOnPatch<Patch1> function_1 = function_values.template get<Patch1>();
        ValuesOnPatch<Patch2> function_2 = function_values.template get<Patch2>();

        auto const idx_range_2d_1 = get_idx_range(function_1);
        auto const idx_range_2d_2 = get_idx_range(function_2);

        auto const idx_range_parell_1 = ddc::remove_dims_of(idx_range_2d_1, idx_range_perp_1);
        auto const idx_range_parell_2 = ddc::remove_dims_of(idx_range_2d_2, idx_range_perp_2);

        auto slice_indexes = get_slice_indexes<
                I,
                Patch1,
                Patch2>(slice_idx_1_value, sorted_idx_ranges_tuple, function_values);
        auto slice_idx_1 = std::get<0>(slice_indexes);
        auto slice_idx_2 = std::get<1>(slice_indexes);

        auto function_slice_1 = function_1[slice_idx_1]; 
        auto function_slice_2 = function_2[slice_idx_2]; 

        double lin_comb_funct = std::get<I>(m_derivatives_calculators)
        .get_function_coefficients(function_slice_1, function_slice_2);
        // =======================================================================================

        if (I == 0) {
            if (LowerBound == ddc::BoundCond::GREVILLE) {
                m_vector->get_values()[I] = fct_coeff;
            } else {
                double coeff_plus = std::get<I>(m_derivatives_calculators).get_coeff_plus();
                auto deriv_min = derivs_min.template get<Patch1>();
                auto idx_range_deriv_1 = ddc::
                        remove_dims_of(get_idx_range(deriv_min), other_idx_range_1d_function_1);

                m_vector->get_values()[I]
                        = fct_coeff
                          + coeff_plus * deriv_min(idx_range_deriv_1.front(), slice_idx_1);
            }



        } else if (1 < I && I < n_inner_interfaces - 1) {
            m_vector->get_values()[I] = fct_coeff;

        } else if (I == n_inner_interfaces - 1) {
            if (UpperBound == ddc::BoundCond::GREVILLE) {
                m_vector->get_values()[I] = fct_coeff;
            } else {
                double coeff_minus = std::get<I>(m_derivatives_calculators).get_coeff_minus();
                auto deriv_max = derivs_max.template get<Patch2>();
                auto idx_range_deriv_2 = ddc::
                        remove_dims_of(get_idx_range(deriv_max), other_idx_range_1d_function_2);

                m_vector->get_values()[I]
                        = fct_coeff
                          + coeff_minus * deriv_max(idx_range_deriv_2.front(), slice_idx_2);
            }
        }
    }


    template <class OGrid1D, std::size_t... I>
    void update_derivatives(
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max,
            std::array<double, n_inner_interfaces> const& derivs_at_interfaces,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            Idx<OGrid1D> const& slice_idx,
            std::integer_sequence<std::size_t, I...>)
    {
        std::tuple sorted_idx_ranges_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);
        using FirstIdxRange = std::tuple_element_t<0, decltype(sorted_idx_ranges_tuple)>;
        using FirstPatch = get_patch_of_idx_range_1d<FirstIdxRange>;
        static_assert(
                (std::is_same_v<OGrid1D, typename FirstPatch::Grid1>)
                || (std::is_same_v<OGrid1D, typename FirstPatch::Grid2>)
                || (std::is_same_v<OGrid1D, ddc::Deriv<typename FirstPatch::Dim1>>)
                || (std::is_same_v<OGrid1D, ddc::Deriv<typename FirstPatch::Dim2>>));

        IdxRange<OGrid1D> idx_range_1d_first_patch(m_idx_ranges.template get<FirstPatch>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (update_derivatives_at_interface<
                 I>(derivs_min, derivs_max, derivs_at_interfaces, function_values, slice_idx_value),
         ...);
    }


    template <std::size_t I>
    void update_derivatives_at_interface(
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max,
            std::array<double, n_inner_interfaces> const& derivs_at_interfaces,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            int& slice_idx_1_value)
    {
        std::tuple sorted_idx_ranges_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);
        static_assert(
                std::tuple_size_v<decltype(sorted_idx_ranges_tuple)> == n_inner_interfaces + 1);

        using Patch1 = get_patch_of_idx_range_1d<
                std::tuple_element_t<I, decltype(sorted_idx_ranges_tuple)>>;
        using Patch2 = get_patch_of_idx_range_1d<
                std::tuple_element_t<I + 1, decltype(sorted_idx_ranges_tuple)>>;

        using interfaces_typeseq =
                typename Connectivity::template find_connections_t<Patch1, Patch2>;
        static_assert(ddc::type_seq_size_v<interfaces_typeseq> == 1);
        using Interface = typename ddc::type_seq_element_t<0, interfaces_typeseq>;

        Extremity const extremity_1 = Interface::Edge1::extremity;
        Extremity const extremity_2 = Interface::Edge2::extremity;

        // Define all the different index ranges
        auto const idx_range_1d_1 = std::get<I>(sorted_idx_ranges_tuple);
        auto const idx_range_1d_2 = std::get<I + 1>(sorted_idx_ranges_tuple);

        auto const function_idx_range_2d_1 = get_idx_range(function_values.template get<Patch1>());
        auto const function_idx_range_2d_2 = get_idx_range(function_values.template get<Patch2>());

        auto const other_idx_range_1d_function_1
                = ddc::remove_dims_of(function_idx_range_2d_1, idx_range_1d_1);
        auto const other_idx_range_1d_function_2
                = ddc::remove_dims_of(function_idx_range_2d_2, idx_range_1d_2);


        auto slice_indexes = get_slice_indexes<
                I,
                Patch1,
                Patch2>(slice_idx_1_value, sorted_idx_ranges_tuple, function_values);
        auto slice_idx_1 = std::get<0>(slice_indexes);
        auto slice_idx_2 = std::get<1>(slice_indexes);

        auto idx_range_deriv_1 = ddc::remove_dims_of(
                get_idx_range(derivs_min.template get<Patch1>()),
                other_idx_range_1d_function_1);
        auto idx_range_deriv_2 = ddc::remove_dims_of(
                get_idx_range(derivs_min.template get<Patch2>()),
                other_idx_range_1d_function_2);


        // Update derivatives
        if (extremity_1 == FRONT) {
            derivs_min.template get<Patch1>()(idx_range_deriv_1.front(), slice_idx_1)
                    = derivs_at_interfaces[I];
        } else {
            derivs_max.template get<Patch1>()(idx_range_deriv_1.front(), slice_idx_1)
                    = derivs_at_interfaces[I];
        }

        if (extremity_2 == FRONT) {
            derivs_min.template get<Patch2>()(idx_range_deriv_2.front(), slice_idx_2)
                    = derivs_at_interfaces[I];
        } else {
            derivs_max.template get<Patch2>()(idx_range_deriv_2.front(), slice_idx_2)
                    = derivs_at_interfaces[I];
        }
    }



    template <std::size_t I, class Patch1, class Patch2, class IdxRangeTuple>
    auto get_slice_indexes(
            int& slice_idx_1_value,
            IdxRangeTuple const& sorted_idx_ranges_tuple,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values)
    {
        using interfaces_typeseq =
                typename Connectivity::template find_connections_t<Patch1, Patch2>;
        static_assert(ddc::type_seq_size_v<interfaces_typeseq> == 1);
        using Interface = typename ddc::type_seq_element_t<0, interfaces_typeseq>;

        // Define index ranges
        typename Patch1::IdxRange12 const idx_range_2d_1 = m_idx_ranges.template get<Patch1>();
        typename Patch2::IdxRange12 const idx_range_2d_2 = m_idx_ranges.template get<Patch2>();

        auto const idx_range_1d_1 = std::get<I>(sorted_idx_ranges_tuple);
        auto const idx_range_1d_2 = std::get<I + 1>(sorted_idx_ranges_tuple);

        auto const other_idx_range_1d_1 = ddc::remove_dims_of(idx_range_2d_1, idx_range_1d_1);
        auto const other_idx_range_1d_2 = ddc::remove_dims_of(idx_range_2d_2, idx_range_1d_2);


        auto const function_idx_range_2d_1 = get_idx_range(function_values.template get<Patch1>());
        auto const function_idx_range_2d_2 = get_idx_range(function_values.template get<Patch2>());

        auto const other_idx_range_1d_function_1
                = ddc::remove_dims_of(function_idx_range_2d_1, idx_range_1d_1);
        auto const other_idx_range_1d_function_2
                = ddc::remove_dims_of(function_idx_range_2d_2, idx_range_1d_2);



        // Get slice indexes
        EdgeTransformation<Interface> index_converter(other_idx_range_1d_1, other_idx_range_1d_2);

        using OIdx1 = typename decltype(other_idx_range_1d_function_1)::discrete_element_type;
        using OIdx2 = typename decltype(other_idx_range_1d_function_2)::discrete_element_type;
        OIdx1 slice_idx_1;
        OIdx2 slice_idx_2;
        if constexpr (std::is_same_v<
                              decltype(other_idx_range_1d_2),
                              decltype(other_idx_range_1d_function_2)>) {
            slice_idx_1 = OIdx1(slice_idx_1_value);
            slice_idx_2 = OIdx2(index_converter(slice_idx_1));
            slice_idx_1_value = (slice_idx_2 - other_idx_range_1d_2.front()).value();

        } else {
            slice_idx_1 = OIdx1(1);
            slice_idx_2 = OIdx2(1);
        }

        std::pair<OIdx1, OIdx2> slice_indexes(slice_idx_1, slice_idx_2);
        return slice_indexes;
    }

    template <
            class Grid_1D,
            class Patch = get_patch_on_grid<Grid_1D>,
            class OIdx = std::conditional_t<
                    std::is_same_v<Grid_1D, typename Patch::Grid1>,
                    typename Patch::Idx2,
                    typename Patch::Idx1>>
    auto const get_function_on_idx_range(
            IdxRange<Grid_1D> const& idx_range_1d,
            OIdx const& idx,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values)
    {
        ValuesOnPatch<Patch> const function_2d = function_values.template get<Patch>();
        return function_2d[idx];
    }
};
