// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "matching_idx_slice.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "single_interface_derivatives_calculator.hpp"
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

    static constexpr std::size_t number_of_interfaces = ddc::type_seq_size_v<interface_collection>;

    // TODO: remove Interfaces with OutsideEdge
    // Remove all the interfaces with an OutsideEdge.
    // Rely on the get_all_interfaces_along_direction_t that orders the interfaces.
    using outer_interface_collection = std::conditional_t<
            PERIODIC,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<
                    ddc::type_seq_element_t<0, interface_collection>,
                    ddc::type_seq_element_t<number_of_interfaces - 1, interface_collection>>>;

    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_collection, outer_interface_collection>;

    using FirstInterface = ddc::type_seq_element_t<0, inner_interface_collection>;
    using LastInterface = ddc::type_seq_element_t<
            ddc::type_seq_size_v<inner_interface_collection> - 1,
            inner_interface_collection>;

    using FirstPatch = std::conditional_t<
            PERIODIC,
            find_patch_t<Grid1D, all_patches>,
            typename FirstInterface::Edge2::associated_patch>;
    using Grid1DSeq = collect_grids_on_dim_t<
            find_patch_t<Grid1D, all_patches>,
            Grid1D,
            interface_collection>;


    template <typename Patch>
    struct GetDerivsOnPatch
    {
        static constexpr bool is_grid_on_dim1 = ddc::in_tags_v<typename Patch::Grid1, Grid1DSeq>;
        using GridPerp
                = std::conditional_t<is_grid_on_dim1, typename Patch::Grid1, typename Patch::Grid2>;
        using GridPar
                = std::conditional_t<is_grid_on_dim1, typename Patch::Grid2, typename Patch::Grid1>;
        using DerivPerp = ddc::Deriv<typename GridPerp::continuous_dimension_type>;
        using DerivPar = ddc::Deriv<typename GridPar::continuous_dimension_type>;

        using IdxRangeDerivPerp = std::conditional_t<
                is_grid_on_dim1,
                IdxRange<DerivPerp, GridPar>,
                IdxRange<GridPar, DerivPerp>>;
        using IdxRangeDerivPar = std::conditional_t<
                is_grid_on_dim1,
                IdxRange<GridPerp, DerivPar>,
                IdxRange<DerivPar, GridPerp>>;

        using type_perp = DField<IdxRangeDerivPerp, typename ExecSpace::memory_space>;
        using type_par = DConstField<IdxRangeDerivPar, typename ExecSpace::memory_space>;
    };

    template <typename Patch>
    using GridPerpOnPatch = typename GetDerivsOnPatch<Patch>::GridPerp;

    template <typename Patch>
    using GridParOnPatch = typename GetDerivsOnPatch<Patch>::GridPar;


    template <typename Patch>
    using DerivsPerpOnPatch = typename GetDerivsOnPatch<Patch>::type_perp;

    template <typename Patch>
    using ConstDerivsParOnPatch = typename GetDerivsOnPatch<Patch>::type_par;

    template <typename Patch>
    using CrossDerivsPerpOnPatch = DField<
            IdxRange<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>>,
            typename ExecSpace::memory_space>;

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



    // Transfrom a tuple type into a TypeSeq.
    template <class Tuple>
    struct TupleToTypeSeq;

    template <class... Elements>
    struct TupleToTypeSeq<std::tuple<Elements...>>
    {
        using type = ddc::detail::TypeSeq<Elements...>;
    };

    template <class Tuple>
    using tuple_to_type_seq_t = typename TupleToTypeSeq<Tuple>::type;


    // ===========================================================================================

    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;

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
                                  gko::stop::ResidualNorm<>::build().with_reduction_factor(1e-6).on(
                                          exec))
                          .with_preconditioner(
                                  PreconditionerJ::build().with_max_block_size(8u).on(exec))
                          .on(exec);
        // Create solver
        m_solver = std::move(solver_factory->generate(m_matrix));

        // Create vectors for rhs and solution
        m_vector = Matrix::create(exec, gko::dim<2>(n_inner_interfaces, 1));
        m_interface_derivatives = Matrix::create(exec, gko::dim<2>(n_inner_interfaces, 1));
    }

    template <class IdxPar>
    void solve(
            IdxPar const& idx_par,
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            MultipatchField<DerivsPerpOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsPerpOnPatch, Patches...> const& derivs_max)
    {
        static_assert(
                (std::is_same_v<IdxPar, typename FirstPatch::Idx1>)
                || (std::is_same_v<IdxPar, typename FirstPatch::Idx2>));
        // Update m_vector.
        set_vector(
                idx_par,
                function_values,
                derivs_min,
                derivs_max,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});

        // Update m_interface_derivatives.
        m_solver->apply(m_vector, m_interface_derivatives);

        // // generate_stencil_matrix(gko::lend(matrix));
        // for (int i = 0; i < n_inner_interfaces; ++i) {
        //     for (int j = 0; j < n_inner_interfaces; ++j) {
        //         std::cout << "matrix = "
        //                   << m_matrix->get_values()[i * (n_inner_interfaces) + j]
        //                   << "   ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;

        // auto values = m_vector->get_values();
        // for (int i = 0; i < n_inner_interfaces; ++i) {
        //     std::cout << "vector = " << values[i] << std::endl;
        // }
        // std::cout << std::endl;

        // auto sol = m_interface_derivatives->get_values();
        // for (int i = 0; i < n_inner_interfaces; ++i) {
        //     std::cout << "sol = " << sol[i] << std::endl;
        // }
        // std::cout << std::endl;

        // Update derivatives.
        update_derivatives(
                derivs_min,
                derivs_max,
                function_values,
                idx_par,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
    }


    void solve(
            MultipatchField<ValuesOnPatch, Patches...> const& function_values,
            MultipatchField<DerivsPerpOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsPerpOnPatch, Patches...> const& derivs_max)
    {
        std::tuple sorted_idx_ranges_perp_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);

        using IdxRangePerpType = tuple_to_type_seq_t<decltype(sorted_idx_ranges_perp_tuple)>;
        constexpr bool is_first_idx_range_perp_on_dim1
                = ddc::in_tags_v<typename FirstPatch::IdxRange1, IdxRangePerpType>;

        using IdxRangePerpFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange1,
                typename FirstPatch::IdxRange2>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange2,
                typename FirstPatch::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first = ddc::remove_dims_of(
                m_idx_ranges.template get<FirstPatch>(),
                std::get<IdxRangePerpFirstType>(sorted_idx_ranges_perp_tuple));

        ddc::for_each(
                idx_range_par_first,
                [&](typename IdxRangeParFirstType::discrete_element_type const& idx_par) {
                    (*this).solve(idx_par, function_values, derivs_min, derivs_max);
                });
    }


    void solve(
            MultipatchField<ConstDerivsParOnPatch, Patches...> const& derivs,
            std::tuple<CrossDerivsPerpOnPatch<Patches>...> const& cross_derivs_min,
            std::tuple<CrossDerivsPerpOnPatch<Patches>...> const& cross_derivs_max)
    {
        std::tuple sorted_idx_ranges_perp_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);

        using IdxRangePerpType = tuple_to_type_seq_t<decltype(sorted_idx_ranges_perp_tuple)>;
        constexpr bool is_first_idx_range_perp_on_dim1
                = ddc::in_tags_v<typename FirstPatch::IdxRange1, IdxRangePerpType>;

        using IdxRangePerpFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange1,
                typename FirstPatch::IdxRange2>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange2,
                typename FirstPatch::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first = ddc::remove_dims_of(
                m_idx_ranges.template get<FirstPatch>(),
                std::get<IdxRangePerpFirstType>(sorted_idx_ranges_perp_tuple));

        ddc::for_each(
                idx_range_par_first,
                [&](typename IdxRangeParFirstType::discrete_element_type const& idx_par) {
                    // Update the m_interface_derivatives vector.
                    set_vector(
                            idx_par,
                            derivs,
                            cross_derivs_min,
                            cross_derivs_max,
                            std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
                    m_solver->apply(m_vector, m_interface_derivatives);

                    update_derivatives(
                            cross_derivs_min,
                            cross_derivs_max,
                            derivs,
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
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        constexpr bool is_same_orientation = std::is_same_v<
                InterfaceI,
                find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>>;

        const double coeff_left
                = is_same_orientation
                          ? std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1()
                          : std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();
        const double coeff_right
                = is_same_orientation
                          ? std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2()
                          : std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
        const double coefs[] = {-coeff_left, 1, -coeff_right};

        for (auto dofs : {-1, 0, 1}) {
            if (0 <= I + dofs && I + dofs < n_inner_interfaces) {
                m_matrix->at(I, I + dofs) = coefs[dofs + 1];
            }
        }

        if constexpr (LowerBound == ddc::BoundCond::PERIODIC && I == 0) {
            m_matrix->at(I, n_inner_interfaces - 1) = coefs[0];
        }
        if constexpr (UpperBound == ddc::BoundCond::PERIODIC && I == n_inner_interfaces - 1) {
            m_matrix->at(I, 0) = coefs[2];
        }
    }

    template <
            class OGrid1D,
            class DerivCollectionType,
            template <typename P>
            class FieldOnPatch,
            std::size_t... I>
    void set_vector(
            Idx<OGrid1D> const& slice_idx,
            MultipatchField<FieldOnPatch, Patches...> const& functions,
            DerivCollectionType const& derivs_min,
            DerivCollectionType const& derivs_max,
            std::integer_sequence<std::size_t, I...>)
    {
        static_assert(
                (std::is_same_v<OGrid1D, typename FirstPatch::Grid1>)
                || (std::is_same_v<OGrid1D, typename FirstPatch::Grid2>)
                || (std::is_same_v<OGrid1D, ddc::Deriv<typename FirstPatch::Dim1>>)
                || (std::is_same_v<OGrid1D, ddc::Deriv<typename FirstPatch::Dim2>>));

        IdxRange<OGrid1D> idx_range_1d_first_patch(m_idx_ranges.template get<FirstPatch>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (set_line_vector<I>(slice_idx_value, functions, derivs_min, derivs_max), ...);
    }


    template <std::size_t I, template <typename P> class FieldOnPatch, class DerivCollectionType>
    void set_line_vector(
            int& slice_idx_1_value,
            MultipatchField<FieldOnPatch, Patches...> const& functions,
            DerivCollectionType const& derivs_min,
            DerivCollectionType const& derivs_max)
    {
        // RMK: only needed for I = 0 and I =  n_inner_interfaces - 1 but computed for each I.
        // auto [deriv_1, deriv_2] = get_interface_derivatives<I>(derivs_min, derivs_min);

        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        constexpr bool is_same_orientation = std::is_same_v<
                InterfaceI,
                find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>>;

        const int sign = is_same_orientation - !is_same_orientation;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        auto [idx_slice_1, idx_slice_2] = get_indexes<I>(slice_idx_1_value, functions);

        FieldOnPatch<Patch1> function_1 = functions.template get<Patch1>();
        FieldOnPatch<Patch2> function_2 = functions.template get<Patch2>();

        auto function_slice_1 = function_1[idx_slice_1];
        auto function_slice_2 = function_2[idx_slice_2];


        // Check that we select the correct derivative calculator.
        // TODO: use a MultipathType over the interfaces?
        using DerivativeCalculator
                = std::decay_t<std::tuple_element_t<I, SingleInterfaceDerivativesCalculatorTuple>>;
        using InterfaceDerivativeCalculator = typename DerivativeCalculator::associated_interface;
        using InterfaceDerivativeCalculatorSymmetrical = Interface<
                typename InterfaceDerivativeCalculator::Edge2,
                typename InterfaceDerivativeCalculator::Edge1,
                InterfaceDerivativeCalculator::orientations_agree>;
        static_assert(
                (is_same_orientation && std::is_same_v<InterfaceI, InterfaceDerivativeCalculator>)
                || (!is_same_orientation
                    && std::is_same_v<InterfaceI, InterfaceDerivativeCalculatorSymmetrical>));

        double lin_comb_funct
                = sign
                  * std::get<I>(m_derivatives_calculators)
                            .get_function_coefficients(function_slice_1, function_slice_2);

        m_vector->get_values()[I] = lin_comb_funct;



        // TODO: Test this part.
        if constexpr (LowerBound == ddc::BoundCond::HERMITE && I == 0) {
            {
                auto [deriv_1, deriv_2]
                        = get_interface_derivatives<I>(derivs_min, derivs_min, Extremity::BACK);
                auto [idx_slice_1, idx_slice_2, idx_first_deriv_1, idx_first_deriv_2]
                        = get_indexes<I>(slice_idx_1_value, functions, deriv_1, deriv_2);

                double coeff_deriv_1
                        = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
                double coeff_deriv_2
                        = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();

                m_vector->get_values()[I]
                        += is_same_orientation
                                   ? coeff_deriv_1 * deriv_1(idx_first_deriv_1, idx_slice_1)
                                   : -coeff_deriv_2 * deriv_2(idx_first_deriv_2, idx_slice_2);
            }
        }
        if constexpr (UpperBound == ddc::BoundCond::HERMITE && I == n_inner_interfaces - 1) {
            auto [deriv_1, deriv_2]
                    = get_interface_derivatives<I>(derivs_min, derivs_min, Extremity::BACK);
            auto [idx_slice_1, idx_slice_2, idx_first_deriv_1, idx_first_deriv_2]
                    = get_indexes<I>(slice_idx_1_value, functions, deriv_1, deriv_2);

            double coeff_deriv_1 = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
            double coeff_deriv_2 = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();

            m_vector->get_values()[I]
                    += is_same_orientation
                               ? coeff_deriv_2 * deriv_2(idx_first_deriv_2, idx_slice_2)
                               : -coeff_deriv_1 * deriv_1(idx_first_deriv_1, idx_slice_1);
        }
    }

    template <
            class OGrid1D,
            class DerivCollectionType,
            template <typename P>
            class FieldOnPatch,
            std::size_t... I>
    void update_derivatives(
            DerivCollectionType const& derivs_min,
            DerivCollectionType const& derivs_max,
            MultipatchField<FieldOnPatch, Patches...> const& function_values,
            Idx<OGrid1D> const& slice_idx,
            std::integer_sequence<std::size_t, I...>)
    {
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


    template <std::size_t I, template <typename P> class FieldOnPatch, class DerivCollectionType>
    void update_derivatives_at_interface(
            DerivCollectionType const& derivs_min,
            DerivCollectionType const& derivs_max,
            MultipatchField<FieldOnPatch, Patches...> const& function_values,
            int& slice_idx_1_value)
    {
        auto [deriv_1, deriv_2]
                = get_interface_derivatives<I>(derivs_min, derivs_max, Extremity::FRONT);

        auto [slice_idx_1, slice_idx_2, idx_first_deriv_1, idx_first_deriv_2]
                = get_indexes<I>(slice_idx_1_value, function_values, deriv_1, deriv_2);

        // Update derivatives
        deriv_1(idx_first_deriv_1, slice_idx_1) = m_interface_derivatives->get_values()[I];
        deriv_2(idx_first_deriv_2, slice_idx_2) = m_interface_derivatives->get_values()[I];
    }



    template <
            std::size_t I,
            class Patch1 = typename ddc::type_seq_element_t<I, inner_interface_collection>::Edge1::
                    associated_patch,
            class Patch2 = typename ddc::type_seq_element_t<I, inner_interface_collection>::Edge2::
                    associated_patch>
    std::tuple<DerivsPerpOnPatch<Patch1>, DerivsPerpOnPatch<Patch2>> get_interface_derivatives(
            MultipatchField<DerivsPerpOnPatch, Patches...> const& derivs_min,
            MultipatchField<DerivsPerpOnPatch, Patches...> const& derivs_max,
            Extremity const extremity)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        constexpr Extremity extermity_1 = InterfaceI::Edge1::extremity;
        constexpr Extremity extermity_2 = InterfaceI::Edge2::extremity;

        DerivsPerpOnPatch<Patch1> deriv_1 = (extermity_1 == extremity)
                                                    ? derivs_min.template get<Patch1>()
                                                    : derivs_max.template get<Patch1>();
        DerivsPerpOnPatch<Patch2> deriv_2 = (extermity_2 == extremity)
                                                    ? derivs_min.template get<Patch2>()
                                                    : derivs_max.template get<Patch2>();
        return std::make_tuple(deriv_1, deriv_2);
    }


    template <
            std::size_t I,
            class Patch1 = typename ddc::type_seq_element_t<I, inner_interface_collection>::Edge1::
                    associated_patch,
            class Patch2 = typename ddc::type_seq_element_t<I, inner_interface_collection>::Edge2::
                    associated_patch>
    std::tuple<CrossDerivsPerpOnPatch<Patch1>, CrossDerivsPerpOnPatch<Patch2>>
    get_interface_derivatives(
            std::tuple<CrossDerivsPerpOnPatch<Patches>...> const& cross_derivs_min,
            std::tuple<CrossDerivsPerpOnPatch<Patches>...> const& cross_derivs_max,
            Extremity const extremity)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        // Hyp: the derivatives in the tuples are well ordered.
        constexpr std::size_t I1 = I;
        constexpr std::size_t I2 = (I + 1) % n_inner_interfaces;

        constexpr Extremity extermity_1 = InterfaceI::Edge1::extremity;
        constexpr Extremity extermity_2 = InterfaceI::Edge2::extremity;

        // TODO: Assert to be sure that the cross-derivatives are well-ordered?
        // How to do it? Dim specialised for patches? message?
        CrossDerivsPerpOnPatch<Patch1> deriv_1 = (extermity_1 == extremity)
                                                         ? std::get<I1>(cross_derivs_min)
                                                         : std::get<I1>(cross_derivs_max);
        CrossDerivsPerpOnPatch<Patch2> deriv_2 = (extermity_2 == extremity)
                                                         ? std::get<I2>(cross_derivs_min)
                                                         : std::get<I2>(cross_derivs_max);
        return std::make_tuple(deriv_1, deriv_2);
    }


    template <std::size_t I, template <typename P> class FieldOnPatch>
    auto get_indexes(
            int& slice_idx_1_value,
            MultipatchField<FieldOnPatch, Patches...> const& functions)
    {
        std::tuple sorted_idx_ranges_perp_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);
        static_assert(
                std::tuple_size_v<
                        decltype(sorted_idx_ranges_perp_tuple)> == n_inner_interfaces + !PERIODIC);

        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        using PerpGrid1 = typename InterfaceI::Edge1::perpendicular_grid;
        using PerpGrid2 = typename InterfaceI::Edge2::perpendicular_grid;

        // Define all the different index ranges
        IdxRange<PerpGrid1> const idx_range_perp_1
                = std::get<IdxRange<PerpGrid1>>(sorted_idx_ranges_perp_tuple);
        IdxRange<PerpGrid2> const idx_range_perp_2
                = std::get<IdxRange<PerpGrid2>>(sorted_idx_ranges_perp_tuple);

        FieldOnPatch<Patch1> function_1 = functions.template get<Patch1>();
        FieldOnPatch<Patch2> function_2 = functions.template get<Patch2>();

        auto const idx_range_fct_2d_1 = get_idx_range(function_1);
        auto const idx_range_fct_2d_2 = get_idx_range(function_2);

        auto const idx_range_fct_parell_1
                = ddc::remove_dims_of(idx_range_fct_2d_1, idx_range_perp_1);
        auto const idx_range_fct_parell_2
                = ddc::remove_dims_of(idx_range_fct_2d_2, idx_range_perp_2);


        // Determine the perpendicular indexes.
        auto slice_indexes = get_slice_indexes<
                I,
                Patch1,
                Patch2>(slice_idx_1_value, sorted_idx_ranges_perp_tuple, functions);
        auto slice_idx_1 = std::get<0>(slice_indexes);
        auto slice_idx_2 = std::get<1>(slice_indexes);

        return std::make_tuple(slice_idx_1, slice_idx_2);
    }

    template <
            std::size_t I,
            template <typename P>
            class FieldOnPatch,
            class DerivFieldOnPatch1,
            class DerivFieldOnPatch2>
    auto get_indexes(
            int& slice_idx_1_value,
            MultipatchField<FieldOnPatch, Patches...> const& functions,
            DerivFieldOnPatch1 const& deriv_1,
            DerivFieldOnPatch2 const& deriv_2)
    {
        std::tuple sorted_idx_ranges_perp_tuple
                = Connectivity::template get_all_idx_ranges_along_direction<Grid1D>(m_idx_ranges);
        static_assert(
                std::tuple_size_v<
                        decltype(sorted_idx_ranges_perp_tuple)> == n_inner_interfaces + !PERIODIC);

        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        using PerpGrid1 = typename InterfaceI::Edge1::perpendicular_grid;
        using PerpGrid2 = typename InterfaceI::Edge2::perpendicular_grid;

        // Define all the different index ranges
        IdxRange<PerpGrid1> const idx_range_perp_1
                = std::get<IdxRange<PerpGrid1>>(sorted_idx_ranges_perp_tuple);
        IdxRange<PerpGrid2> const idx_range_perp_2
                = std::get<IdxRange<PerpGrid2>>(sorted_idx_ranges_perp_tuple);

        FieldOnPatch<Patch1> function_1 = functions.template get<Patch1>();
        FieldOnPatch<Patch2> function_2 = functions.template get<Patch2>();

        auto const idx_range_fct_2d_1 = get_idx_range(function_1);
        auto const idx_range_fct_2d_2 = get_idx_range(function_2);

        auto const idx_range_fct_parell_1
                = ddc::remove_dims_of(idx_range_fct_2d_1, idx_range_perp_1);
        auto const idx_range_fct_parell_2
                = ddc::remove_dims_of(idx_range_fct_2d_2, idx_range_perp_2);


        // Determine the perpendicular indexes.
        auto slice_indexes = get_slice_indexes<
                I,
                Patch1,
                Patch2>(slice_idx_1_value, sorted_idx_ranges_perp_tuple, functions);
        auto slice_idx_1 = std::get<0>(slice_indexes);
        auto slice_idx_2 = std::get<1>(slice_indexes);

        // Determine the parallel indexes.
        Idx<ddc::Deriv<typename PerpGrid1::continuous_dimension_type>> idx_first_deriv_1
                = ddc::remove_dims_of(get_idx_range(deriv_1), idx_range_fct_parell_1).front();
        Idx<ddc::Deriv<typename PerpGrid2::continuous_dimension_type>> idx_first_deriv_2
                = ddc::remove_dims_of(get_idx_range(deriv_2), idx_range_fct_parell_2).front();

        return std::make_tuple(slice_idx_1, slice_idx_2, idx_first_deriv_1, idx_first_deriv_2);
    }


    // Get the equivalent index slice on the two patches.
    template <
            std::size_t I,
            class Patch1,
            class Patch2,
            class IdxRangeTuple,
            template <typename P>
            class FieldOnPatch>
    auto get_slice_indexes(
            int& slice_idx_1_value,
            IdxRangeTuple const& sorted_idx_ranges_tuple,
            MultipatchField<FieldOnPatch, Patches...> const& function_values)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        using PerpGrid1 = typename InterfaceI::Edge1::perpendicular_grid;
        using PerpGrid2 = typename InterfaceI::Edge2::perpendicular_grid;

        using ParallGrid1 = typename InterfaceI::Edge1::parallel_grid;
        using ParallGrid2 = typename InterfaceI::Edge2::parallel_grid;

        // Define all the different index ranges
        IdxRange<PerpGrid1> const idx_range_perp_1
                = std::get<IdxRange<PerpGrid1>>(sorted_idx_ranges_tuple);
        IdxRange<PerpGrid2> const idx_range_perp_2
                = std::get<IdxRange<PerpGrid2>>(sorted_idx_ranges_tuple);

        typename Patch1::IdxRange12 const idx_range_2d_1 = m_idx_ranges.template get<Patch1>();
        typename Patch2::IdxRange12 const idx_range_2d_2 = m_idx_ranges.template get<Patch2>();

        IdxRange<ParallGrid1> const idx_range_parell_1
                = ddc::remove_dims_of(idx_range_2d_1, idx_range_perp_1);
        IdxRange<ParallGrid2> const idx_range_parell_2
                = ddc::remove_dims_of(idx_range_2d_2, idx_range_perp_2);


        FieldOnPatch<Patch1> function_1 = function_values.template get<Patch1>();
        FieldOnPatch<Patch2> function_2 = function_values.template get<Patch2>();

        auto const idx_range_fct_2d_1 = get_idx_range(function_1);
        auto const idx_range_fct_2d_2 = get_idx_range(function_2);

        auto const idx_range_fct_parall_1
                = ddc::remove_dims_of(idx_range_fct_2d_1, idx_range_perp_1);
        auto const idx_range_fct_parall_2
                = ddc::remove_dims_of(idx_range_fct_2d_2, idx_range_perp_2);

        // Get slice indexes
        EdgeTransformation<InterfaceI> index_converter(idx_range_parell_1, idx_range_parell_2);

        using OIdx1 = typename decltype(idx_range_fct_parall_1)::discrete_element_type;
        using OIdx2 = typename decltype(idx_range_fct_parall_2)::discrete_element_type;
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
