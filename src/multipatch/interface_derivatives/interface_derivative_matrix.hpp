// SPDX-License-Identifier: MIT

#pragma once

// #include <Eigen>

#include <ddc/ddc.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "matching_idx_slice.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "single_interface_derivatives_calculator.hpp"
// #include "geometry_descriptor.hpp"
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
 *          - update derivatives. How to assert values?
 *
 *  - Remarks:
 *      - Should we restrain to the conforming case first?
 *      - And then the non-conforming case in another operator?
 *      - Treat only "one block line/direction"?
 *      - And build another operator to manage the whole geometry?
 *      - Including treatment for T-joints?
 *      - Deal with periodic case here or in another operator?
 */

template <
        class Connectivity,
        class Grid1D,
        template <typename P>
        typename ValuesOnPatch,
        template <typename P>
        typename DerivsOnPatch,
        bool PERIODIC,
        ddc::BoundCond LowerBound = ddc::BoundCond::HERMITE,
        ddc::BoundCond UpperBound = ddc::BoundCond::HERMITE,
        class ExecSpace = Kokkos::DefaultHostExecutionSpace,
        class... Patches>
class InterfaceDerivativeMatrix
{
    using all_interface_collection = typename Connectivity::interface_collection; // TypeSeq
    using interface_collection =
            typename Connectivity::get_all_interfaces_along_direction_t<Grid1D>;

    using all_patches = typename Connectivity::all_patches; // TypeSeq

    // using grid_collection = collect_grids_on_dim_t<Patch1, Grid1D, interface_collection>;

    // TODO: remove Interfaces with OutsideEdge
    // using inner_interface_collection = ddc::detail::TypeSeq<
    //         ddc::type_seq_element_t<8, interface_collection>,
    //         ddc::type_seq_element_t<9, interface_collection>>;

    static constexpr std::size_t number_of_interfaces = ddc::type_seq_size_v<interface_collection>;

    // Remove all the interfaces with an OutsideEdge.
    using outer_interface_collection = std::conditional_t<
            PERIODIC,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<
                    ddc::type_seq_element_t<0, interface_collection>,
                    ddc::type_seq_element_t<number_of_interfaces - 1, interface_collection>>>;

    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_collection, outer_interface_collection>;

    // using inner_inner_interface_collection = ddc::type_seq_remove_t<
    //         ddc::type_seq_element_t<1, interface_collection>,
    //         ddc::type_seq_element_t<number_of_interfaces - 2, interface_collection>,
    //         inner_interface_collection>;


    // HELPFUL ALIASES ===========================================================================
    template <typename Grid>
    using get_patch_on_grid = find_patch_t<Grid, all_patches>;


    // Get the grid tag from a tag templated on the grid.
    template <typename T>
    struct get_grid;

    template <template <typename P> typename T, typename Grid>
    struct get_grid<T<Grid>>
    {
        using type = Grid;
    };

    template <typename T>
    using get_grid_t = typename get_grid<T>::type;

    // Get the Patch defined with the given IdxRange tag.
    template <typename IdxRange1D>
    using get_patch_of_idx_range_1d = get_patch_on_grid<get_grid_t<IdxRange1D>>;

    // Get a tuple of tags templated on the Interface with the order of the
    // given Interface collection.
    template <template <typename P> typename T, class InterfacesTypeSeq>
    struct get_tuple_on_interfaces;

    template <template <typename P> typename T, class InterfacesTypeSeq>
    using get_tuple_on_interfaces_t = typename get_tuple_on_interfaces<T, InterfacesTypeSeq>::type;

    template <template <typename P> typename T, class... Interfaces>
    struct get_tuple_on_interfaces<T, ddc::detail::TypeSeq<Interfaces...>>
    {
        using type = std::tuple<T<Interfaces> const&...>;
    };


    // Get the type of the interface given to define the Connectivity class.
    template <typename CurrentInterface>
    using get_equivalent_interface_t = find_associated_interface_t<
            typename CurrentInterface::Edge1,
            all_interface_collection>;


    // Get a tuple of SingleInterfaceDerivativesCalculator from an Interface collection.
    template <bool is_periodic, typename InterfaceCollection>
    struct get_tuple_deriv_calculator;

    template <class... Interfaces>
    struct get_tuple_deriv_calculator<true, ddc::detail::TypeSeq<Interfaces...>>
    {
        using InterfaceCollection = ddc::detail::TypeSeq<Interfaces...>;
        using type = std::tuple<SingleInterfaceDerivativesCalculator<
                get_equivalent_interface_t<Interfaces>> const&...>;
    };

    template <class... Interfaces>
    struct get_tuple_deriv_calculator<false, ddc::detail::TypeSeq<Interfaces...>>
    {
        using InterfaceCollection = ddc::detail::TypeSeq<Interfaces...>;

        using FirstInterface = ddc::type_seq_element_t<0, InterfaceCollection>;
        using LastInterface = ddc::type_seq_element_t<
                ddc::type_seq_size_v<InterfaceCollection> - 1,
                InterfaceCollection>;

        using EquivalentFirstInterface = get_equivalent_interface_t<FirstInterface>;
        using EquivalentLastInterface = get_equivalent_interface_t<LastInterface>;

        static constexpr bool is_first_interface_same_orientation
                = std::is_same_v<FirstInterface, EquivalentFirstInterface>;
        static constexpr bool is_last_interface_same_orientation
                = std::is_same_v<LastInterface, EquivalentLastInterface>;

        using inner_inner_interface_collection = ddc::type_seq_remove_t<
                InterfaceCollection,
                ddc::detail::TypeSeq<FirstInterface, LastInterface>>;


        using inner_deriv_calculators =
                typename get_tuple_deriv_calculator<true, inner_inner_interface_collection>::type;

        using first_deriv_calculator = std::conditional_t<
                is_first_interface_same_orientation,
                std::tuple<SingleInterfaceDerivativesCalculator<
                        FirstInterface,
                        LowerBound,
                        ddc::BoundCond::HERMITE> const&>,
                std::tuple<SingleInterfaceDerivativesCalculator<
                        EquivalentFirstInterface,
                        ddc::BoundCond::HERMITE,
                        LowerBound> const&>>;

        using last_deriv_calculator = std::conditional_t<
                is_last_interface_same_orientation,
                std::tuple<SingleInterfaceDerivativesCalculator<
                        LastInterface,
                        ddc::BoundCond::HERMITE,
                        UpperBound> const&>,
                std::tuple<SingleInterfaceDerivativesCalculator<
                        EquivalentLastInterface,
                        UpperBound,
                        ddc::BoundCond::HERMITE> const&>>;

        using type = decltype(std::tuple_cat(
                std::declval<first_deriv_calculator>(),
                std::declval<inner_deriv_calculators>(),
                std::declval<last_deriv_calculator>()));
    };

    template <typename InterfaceTypSeq>
    using get_tuple_deriv_calculator_t =
            typename get_tuple_deriv_calculator<PERIODIC, InterfaceTypSeq>::type;



    // ===========================================================================================

    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;
    // static constexpr std::size_t n_values = n_inner_interfaces * 3 - 2; // TODO: if periodic, we do remove 2.

    using Matrix = gko::matrix::Dense<double>;



private:
    // Matrix m_matrix;
    std::shared_ptr<Matrix> m_matrix;
    std::shared_ptr<Matrix> m_vector;
    std::shared_ptr<Matrix> m_interface_derivatives;


    // Use a conjugate gradient (CG) solver.
    using SolverCG = gko::solver::Cg<>;
    // Use a Jacobi preconditioner.
    using PreconditionerJ = gko::preconditioner::Jacobi<>;
    std::unique_ptr<SolverCG> m_solver;

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


    // WARNING, BE ABLE TO DEAL WITH DIFFERENT ORIENTATION INTERFACES.
    using SingleInterfaceDerivativesCalculatorTuple
            = get_tuple_deriv_calculator_t<inner_interface_collection>;
    SingleInterfaceDerivativesCalculatorTuple const& m_derivatives_calculators;

public:
    ~InterfaceDerivativeMatrix() {}

    InterfaceDerivativeMatrix(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            SingleInterfaceDerivativesCalculatorTuple const derivatives_calculators)
        : m_idx_ranges(idx_ranges)
        , m_derivatives_calculators(derivatives_calculators)
    {
        const auto exec = gko::ReferenceExecutor::create();
        m_matrix = gko::share(Matrix::create(exec, gko::dim<2>(n_inner_interfaces)));

        set_matrix(std::make_integer_sequence<std::size_t, n_inner_interfaces> {});

        // Define a solver.
        auto solver_factory
                = SolverCG::build()
                          .with_criteria(
                                  gko::stop::Iteration::build()
                                          .with_max_iters(n_inner_interfaces)
                                          .on(exec),
                                  gko::stop::ResidualNormReduction<>::build()
                                          .with_reduction_factor(1e-6)
                                          .on(exec))
                          .with_preconditioner(
                                  PreconditionerJ::build().with_max_block_size(8u).on(exec))
                          .on(exec);
        // Create solver
        m_solver = std::move(solver_factory->generate(m_matrix));
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
        auto idx_range_par_0 = ddc::remove_dims_of(
                m_idx_ranges.template get<FirstPatch>(),
                std::get<0>(sorted_idx_ranges_tuple));

        ddc::for_each(idx_range_par_0, [&](auto const& idx_par) {
            // Update the m_interface_derivatives vector.
            set_vector(
                    idx_par,
                    function_values,
                    derivs_min,
                    derivs_max,
                    std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
            m_solver->apply(m_vector, m_interface_derivatives);


            update_derivatives(
                    derivs_min,
                    derivs_max,
                    m_interface_derivatives,
                    function_values,
                    idx_par,
                    std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
        });
    }


private:
    template <std::size_t... I>
    void set_matrix(std::integer_sequence<std::size_t, I...>)
    {
        (set_line_matrix<I>(), ...);
    }

    // Fill in the line of the matrix corresponding to the given interface.
    template <std::size_t I>
    void set_line_matrix()
    {
        const double coeff_right = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
        const double coeff_left = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();
        const double coefs[] = {-coeff_left, 1, -coeff_right};

        for (auto dofs : {-1, 0, 1}) {
            if (0 <= I + dofs && I + dofs < n_inner_interfaces) {
                m_matrix->at(I, I + dofs) = coefs[dofs + 1];
            }
        }
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

        // TODO: Are the Interface sorted? Is it possible to have them sorted?
        using Interface = ddc::type_seq_element_t<I, interface_collection>;
        // using Patch1 = typename  Interface::Edge1::associated_patch;
        // using Patch2 = typename  Interface::Edge2::associated_patch;
        const Extremity extermity_1 = Interface::Edge1::extermity;
        const Extremity extermity_2 = Interface::Edge2::extermity;

        using PerpGrid1 = typename Interface::Edge1::perpendicular_grid;
        using PerpGrid2 = typename Interface::Edge2::perpendicular_grid;

        using ParallGrid1 = typename Interface::Edge1::parallel_grid;
        using ParallGrid2 = typename Interface::Edge2::parallel_grid;


        IdxRange<PerpGrid1> const idx_range_perp_1 = std::get<I>(sorted_idx_ranges_tuple);
        IdxRange<PerpGrid2> const idx_range_perp_2 = std::get<I + 1>(sorted_idx_ranges_tuple);

        ValuesOnPatch<Patch1> function_1 = function_values.template get<Patch1>();
        ValuesOnPatch<Patch2> function_2 = function_values.template get<Patch2>();

        // Remark: not sure that it is:
        // typename Patch1::IdxRange12
        // typename Patch2::IdxRange12
        auto const idx_range_fct_2d_1 = get_idx_range(function_1);
        auto const idx_range_fct_2d_2 = get_idx_range(function_2);

        // Same:
        // IdxRange<ParallGrid1>
        // IdxRange<ParallGrid2>
        auto const idx_range_fct_parell_1
                = ddc::remove_dims_of(idx_range_fct_2d_1, idx_range_perp_1);
        auto const idx_range_fct_parell_2
                = ddc::remove_dims_of(idx_range_fct_2d_2, idx_range_perp_2);

        auto slice_indexes = get_slice_indexes<
                I,
                Patch1,
                Patch2>(slice_idx_1_value, sorted_idx_ranges_tuple, function_values);
        Idx<ParallGrid1> idx_slice_1 = std::get<0>(slice_indexes);
        Idx<ParallGrid2> idx_slice_2 = std::get<1>(slice_indexes);

        auto function_slice_1 = function_1[idx_slice_1];
        auto function_slice_2 = function_2[idx_slice_2];

        double lin_comb_funct
                = std::get<I>(m_derivatives_calculators)
                          .get_function_coefficients(function_slice_1, function_slice_2);


        if (I == 0) {
            // v_I = c + a*deriv_1
            if constexpr (LowerBound == ddc::BoundCond::GREVILLE) {
                m_vector->get_values()[I] = lin_comb_funct;
            } else {
                double coeff_deriv_1
                        = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
                DerivsOnPatch<Patch1> deriv_1 = (extermity_1 == Extremity::BACK)
                                                        ? derivs_min.template get<Patch1>()
                                                        : derivs_max.template get<Patch1>();
                // Idx<ddc::Deriv<>>?
                Idx<ddc::Deriv<typename PerpGrid1::continuous_dimension_type>> idx_first_deriv_1
                        = ddc::remove_dims_of(get_idx_range(deriv_1), idx_range_fct_parell_1)
                                  .front();

                m_vector->get_values()[I]
                        = lin_comb_funct + coeff_deriv_1 * deriv_1(idx_first_deriv_1, idx_slice_1);
            }
        } else if (1 < I && I < n_inner_interfaces - 1) {
            // v_I = c
            m_vector->get_values()[I] = lin_comb_funct;
        } else if (I == n_inner_interfaces - 1) {
            // v_I = c + b*deriv_2
            if constexpr (UpperBound == ddc::BoundCond::GREVILLE) {
                m_vector->get_values()[I] = lin_comb_funct;
            } else {
                double coeff_deriv_2
                        = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();
                DerivsOnPatch<Patch2> deriv_2 = (extermity_2 == Extremity::BACK)
                                                        ? derivs_min.template get<Patch2>()
                                                        : derivs_max.template get<Patch2>();
                Idx<ddc::Deriv<typename PerpGrid2::continuous_dimension_type>> idx_first_deriv_2
                        = ddc::remove_dims_of(get_idx_range(deriv_2), idx_range_fct_parell_2)
                                  .front();

                m_vector->get_values()[I]
                        = lin_comb_funct + coeff_deriv_2 * deriv_2(idx_first_deriv_2, idx_slice_2);
            }
        }
    }


    template <class OGrid1D, std::size_t... I>
    void update_derivatives(
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max,
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
                 I>(derivs_min, derivs_max, function_values, slice_idx_value),
         ...);
    }


    template <std::size_t I>
    void update_derivatives_at_interface(
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsOnPatch, Patches...> const& derivs_max,
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

        using PerpGrid1 = typename Interface::Edge1::perpendicular_grid;
        using PerpGrid2 = typename Interface::Edge2::perpendicular_grid;

        // using ParallGrid1 = typename Interface::Edge1::parallel_grid;
        // using ParallGrid2 = typename Interface::Edge2::parallel_grid;

        // Define all the different index ranges
        IdxRange<PerpGrid1> const idx_range_perp_1 = std::get<I>(sorted_idx_ranges_tuple);
        IdxRange<PerpGrid2> const idx_range_perp_2 = std::get<I + 1>(sorted_idx_ranges_tuple);

        ValuesOnPatch<Patch1> function_1 = function_values.template get<Patch1>();
        ValuesOnPatch<Patch2> function_2 = function_values.template get<Patch2>();

        // Remark: not sure that it is:
        // typename Patch1::IdxRange12
        // typename Patch2::IdxRange12
        auto const idx_range_fct_2d_1 = get_idx_range(function_1);
        auto const idx_range_fct_2d_2 = get_idx_range(function_2);

        // Same:
        // IdxRange<ParallGrid1>
        // IdxRange<ParallGrid2>
        auto const idx_range_fct_parell_1
                = ddc::remove_dims_of(idx_range_fct_2d_1, idx_range_perp_1);
        auto const idx_range_fct_parell_2
                = ddc::remove_dims_of(idx_range_fct_2d_2, idx_range_perp_2);


        auto slice_indexes = get_slice_indexes<
                I,
                Patch1,
                Patch2>(slice_idx_1_value, sorted_idx_ranges_tuple, function_values);
        auto slice_idx_1 = std::get<0>(slice_indexes);
        auto slice_idx_2 = std::get<1>(slice_indexes);

        DerivsOnPatch<Patch1> deriv_1 = (extremity_1 == Extremity::BACK)
                                                ? derivs_min.template get<Patch1>()
                                                : derivs_max.template get<Patch1>();
        DerivsOnPatch<Patch2> deriv_2 = (extremity_2 == Extremity::BACK)
                                                ? derivs_min.template get<Patch2>()
                                                : derivs_max.template get<Patch2>();

        Idx<ddc::Deriv<typename PerpGrid1::continuous_dimension_type>> idx_first_deriv_1
                = ddc::remove_dims_of(get_idx_range(deriv_1), idx_range_fct_parell_1).front();
        Idx<ddc::Deriv<typename PerpGrid2::continuous_dimension_type>> idx_first_deriv_2
                = ddc::remove_dims_of(get_idx_range(deriv_2), idx_range_fct_parell_2).front();


        // Update derivatives
        deriv_1(idx_first_deriv_1, slice_idx_1) = m_interface_derivatives->get_values()[I];
        deriv_2(idx_first_deriv_2, slice_idx_2) = m_interface_derivatives->get_values()[I];
    }


    // Get the equivalent index slice on the two patches.
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


        using PerpGrid1 = typename Interface::Edge1::perpendicular_grid;
        using PerpGrid2 = typename Interface::Edge2::perpendicular_grid;

        using ParallGrid1 = typename Interface::Edge1::parallel_grid;
        using ParallGrid2 = typename Interface::Edge2::parallel_grid;

        // Define all the different index ranges
        IdxRange<PerpGrid1> const idx_range_perp_1 = std::get<I>(sorted_idx_ranges_tuple);
        IdxRange<PerpGrid2> const idx_range_perp_2 = std::get<I + 1>(sorted_idx_ranges_tuple);

        typename Patch1::IdxRange12 const idx_range_2d_1 = m_idx_ranges.template get<Patch1>();
        typename Patch2::IdxRange12 const idx_range_2d_2 = m_idx_ranges.template get<Patch2>();

        IdxRange<ParallGrid1> const idx_range_parell_1
                = ddc::remove_dims_of(idx_range_2d_1, idx_range_perp_1);
        IdxRange<ParallGrid2> const idx_range_parell_2
                = ddc::remove_dims_of(idx_range_2d_2, idx_range_perp_2);


        ValuesOnPatch<Patch1> function_1 = function_values.template get<Patch1>();
        ValuesOnPatch<Patch2> function_2 = function_values.template get<Patch2>();

        auto const idx_range_fct_2d_1 = get_idx_range(function_1);
        auto const idx_range_fct_2d_2 = get_idx_range(function_2);

        auto const idx_range_fct_parall_1
                = ddc::remove_dims_of(idx_range_fct_2d_1, idx_range_perp_1);
        auto const idx_range_fct_parall_2
                = ddc::remove_dims_of(idx_range_fct_2d_2, idx_range_perp_2);

        // Get slice indexes
        EdgeTransformation<Interface> index_converter(idx_range_parell_1, idx_range_parell_2);

        using OIdx1 = typename decltype(idx_range_fct_parall_1)::discrete_element_type;
        using OIdx2 = typename decltype(idx_range_fct_parall_1)::discrete_element_type;
        OIdx1 slice_idx_1;
        OIdx2 slice_idx_2;
        /*
            If function is defined on the same index range than the ones given,
            then it corresponds to the values of the function we are working on.
            Otherwise, it corresponds to the first derivatives of the function
            we are working. The operator is applied to compute the cross-derivatives.
        */
        if constexpr (std::is_same_v<
                              decltype(idx_range_parell_1),
                              decltype(idx_range_fct_parall_1)>) {
            slice_idx_1 = OIdx1(slice_idx_1_value);
            slice_idx_2 = OIdx2(index_converter(slice_idx_1));
            slice_idx_1_value = (slice_idx_2 - idx_range_parell_2.front()).value();

        } else {
            // Stay index for the first derivatives.
            slice_idx_1 = OIdx1(1);
            slice_idx_2 = OIdx2(1);
        }

        std::pair<OIdx1, OIdx2> slice_indexes(slice_idx_1, slice_idx_2);
        return slice_indexes;
    }
};
