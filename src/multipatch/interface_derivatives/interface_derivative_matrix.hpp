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
#include "view.hpp"

/**
 * TODO:
 *  - How does it work?
 *      - Constructor
 *          - Take a list of "parallel" interfaces.
 *          - For each "conforming" lines, get the coefficients a and b.
 *  - Remarks:
 *      - Should we restrain to the conforming case first?
 *      - And then the non-conforming case in another operator?
 *      - Treat only "one block line/direction"?
 *      - And build another operator to manage the whole geometry?
 *      - Including treatment for T-joints?
 *      - Deal with periodic case here or in another operator?
 */


/**
  * @brief Class to compute the interface derivatives along a given direction 
  * on all the interfaces.  
  * 
  * This operator is implemented for conforming equivalent global meshes. 
  * 
  * During the instantiation of the class, the matrix (I-M) is computed are stored in the class. 
  * When we call the operator .solve(), 
  *     - InterfaceDerivativeMatrix computes the vector C from the given values. 
  *     - It inverses the matrix system: (I-M)S = C. 
  *     - The interface derivatives are stored in the vector S. 
  *     - It updates the given derivatives field with the values of the vector S. 
  * 
  * 
  * @tparam Connectivity A MultipatchConnectivity class describing all the patch connections.
  * @tparam Grid1D A given direction.
  * @tparam ValuesOnPatch Type of the const Field  where the function values are defined. 
  *         It is templated on the patch. 
  * @tparam LowerBound Lower/left boundary condition of the first local spline along the given direction. 
  * @tparam UpperBound Upper/right boundary condition of the last local spline along the given direction. 
  * @tparam ExecSpace Execution space.
  * @tparam Patches List of patches containing all the involved patches in the given direction. 
  */
template <
        class Connectivity,
        class Grid1D,
        class SingleInterfaceDerivativesCalculatorCollectionType,
        ddc::BoundCond LowerBound = ddc::BoundCond::HERMITE,
        ddc::BoundCond UpperBound = ddc::BoundCond::HERMITE,
        class ExecSpace = Kokkos::DefaultHostExecutionSpace,
        class... Patches>
class InterfaceDerivativeMatrix
{
    // All the interfaces given as input to the MultipatchConnectivity class.
    // We expect the parameters defined on these interfaces.
    using all_interface_collection = typename Connectivity::interface_collection;
    // All the sorted interfaces with the correct orientation.
    using interface_collection =
            typename Connectivity::get_all_interfaces_along_direction_t<Grid1D>;

    // All the patches of the geometry.
    using all_patches = typename Connectivity::all_patches;

    static constexpr std::size_t number_of_interfaces = ddc::type_seq_size_v<interface_collection>;

    static_assert(
            ((LowerBound == ddc::BoundCond::PERIODIC) == (UpperBound == ddc::BoundCond::PERIODIC)),
            "If one boundary is periodic, the other boundary should be too.");

    static constexpr bool is_periodic = (LowerBound == ddc::BoundCond::PERIODIC);

    // TODO: remove Interfaces with OutsideEdge
    // Remove all the interfaces with an OutsideEdge.
    // Rely on the get_all_interfaces_along_direction_t that orders the interfaces.
    using outer_interface_collection = std::conditional_t<
            is_periodic,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<
                    ddc::type_seq_element_t<0, interface_collection>,
                    ddc::type_seq_element_t<number_of_interfaces - 1, interface_collection>>>;

    // All the interfaces in the Grid1D direction without the interfaces with OutsideEdge.
    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_collection, outer_interface_collection>;

    // Number of inner interfaces. It fixes the number of element in the matrix and vectors.
    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;

    using FirstInterface = ddc::type_seq_element_t<0, inner_interface_collection>;
    using LastInterface = ddc::type_seq_element_t<
            ddc::type_seq_size_v<inner_interface_collection> - 1,
            inner_interface_collection>;

    // First patch of the sorted interfaces.
    using FirstPatch = std::conditional_t<
            is_periodic,
            find_patch_t<Grid1D, all_patches>,
            typename FirstInterface::Edge2::associated_patch>;

    // Sequence of 1D grids in the direction of Grid1D over the patches.
    using Grid1DSeq = collect_grids_on_dim_t<
            find_patch_t<Grid1D, all_patches>,
            Grid1D,
            interface_collection>;

    static_assert(
            ddc::type_seq_size_v<Grid1DSeq> == n_inner_interfaces + !is_periodic,
            "The number of 1D grids and inner interfaces should match.");

    static_assert(
            ddc::type_seq_size_v<Grid1DSeq> <= sizeof...(Patches),
            "The number of given patches is insufficient.");

    static_assert(
            ddc::type_seq_contains_v<
                    Grid1DSeq,
                    ddc::type_seq_merge_t<
                            ddc::detail::TypeSeq<typename Patches::Grid1...>,
                            ddc::detail::TypeSeq<typename Patches::Grid2...>>>,
            "A patch is missing in the given patch list. The list has to contain all the patches "
            "in the given Grid1D direction.");

    // Struct to define the type of derivatives and cross-derivatives for each patch in the direction.
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

    // Field for the first derivatives along the perpendicular direction of the interface.
    template <typename Patch>
    using DerivsPerpOnPatch = typename GetDerivsOnPatch<Patch>::type_perp;

    // template <class Patch>
    // using DerivFieldOnPatch = DerivField<
    //         double,
    //         IdxRange<
    //                 ddc::Deriv<typename Patch::Dim1>,
    //                 typename Patch::Grid1,
    //                 ddc::Deriv<typename Patch::Dim2>,
    //                 typename Patch::Grid2>>;

    // template <class Patch>
    // using DerivFieldOnPatch_host = host_t<DerivFieldOnPatch<Patch>>;


    // template <class Patch>
    // using IdxRange1SliceOnPatch = IdxRangeSlice<typename Patch::Grid1>;

    // template <class Patch>
    // using IdxRange2SliceOnPatch = IdxRangeSlice<typename Patch::Grid2>;

    struct eval_deriv
    {
    };
    struct eval_cross_deriv
    {
    };


    // HELPFUL ALIASES ===========================================================================
    // Get the type of the interface given to define the Connectivity class.
    template <typename CurrentInterface>
    using get_equivalent_interface_t = find_associated_interface_t<
            typename CurrentInterface::Edge1,
            all_interface_collection>;

    // Get the type of the interface sorted along the Grid1D direction.
    // template <typename CurrentInterface>
    // using get_sorted_interface_t
    //         = find_associated_interface_t<typename CurrentInterface::Edge1, interface_collection>;

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

    template <typename InterfaceTypeSeq>
    using get_tuple_deriv_calculator_t =
            typename get_tuple_deriv_calculator<is_periodic, InterfaceTypeSeq>::type;

    // ===========================================================================================

    using Matrix = gko::matrix::Dense<double>;

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using first_patch = FirstPatch;

private:
    std::shared_ptr<Matrix> m_matrix;
    std::shared_ptr<Matrix> m_vector;
    std::shared_ptr<Matrix> m_interface_derivatives;


    // Use a conjugate gradient (CG) solver.
    using SolverCG = gko::solver::Cg<>;
    // Use a Jacobi preconditioner.
    using PreconditionerJ = gko::preconditioner::Jacobi<>;
    std::unique_ptr<SolverCG> m_solver;

    MultipatchType<IdxRangeOnPatch, Patches...> const& m_idx_ranges;

    // TODO: WARNING, BE ABLE TO DEAL WITH DIFFERENT ORIENTATION/ORDER INTERFACES.
    using SingleInterfaceDerivativesCalculatorTuple
            = get_tuple_deriv_calculator_t<inner_interface_collection>;
    // SingleInterfaceDerivativesCalculatorTuple const& m_derivatives_calculators;

    SingleInterfaceDerivativesCalculatorCollectionType const& m_derivatives_calculators;

public:
    ~InterfaceDerivativeMatrix() {}

    /**
     * @brief Instantiate InterfaceDerivativeMatrix. 
     * 
     * It creates a dense Gingko matrix, a vector for the rhs and a vector for the solution. 
     * The matrix is initialised with the derivatives calculator collection given in input. 
     * 
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param derivatives_calculators SingleInterfaceDerivativesCalculatorCollection containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     * @param reduction_factor Reduction factor for a Gingko dense matrix with a conjugate gradient solver 
     *          and a Jacobi preconditioner. 
     * @warning The interface derivative calculators have to be ordered in the tuple. 
     *      The order has to be the same as the interfaces in the Grid1D direction. 
     */
    InterfaceDerivativeMatrix(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            SingleInterfaceDerivativesCalculatorCollectionType const derivatives_calculators,
            double const& reduction_factor = 1e-6)
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
                                  gko::stop::ResidualNorm<>::build()
                                          .with_reduction_factor(reduction_factor)
                                          .on(exec))
                          .with_preconditioner(
                                  PreconditionerJ::build().with_max_block_size(8u).on(exec))
                          .on(exec);
        // Create solver
        m_solver = std::move(solver_factory->generate(m_matrix));

        // Create vectors for rhs and solution
        m_vector = Matrix::create(exec, gko::dim<2>(n_inner_interfaces, 1));
        m_interface_derivatives = Matrix::create(exec, gko::dim<2>(n_inner_interfaces, 1));
    }



    /**
     * @brief Solve the matrix system MS = C to determine all the interface derivatives in the Grid1D
     * direction at a given index.
     * 
     * It uses the function values to compute the first-derivatives perpendicular to the interfaces. 
     *
     * @tparam IdxPar Type for the given index. It can be on the first or the second grid of the first patch.
     *
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @param[in] idx_par Index of the line in the Grid1D direction where we compute all the interface derivatives.
     */
    template <class IdxPar>
    void solve_deriv( // Should be useful for non-conforming case (later).
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        solve<eval_deriv>(functions_and_derivs, idx_par);
    }

    /**
     * @brief Solve the matrix system MS = C to determine all the interface derivatives in the Grid1D
     * direction at all the indices of the first patch.
     * 
     * It uses the function values to compute the first-derivatives perpendicular to the interfaces. 
     * 
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives 
     *              we want to compute.
     */
    void solve_deriv(MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        constexpr bool is_first_idx_range_perp_on_dim1
                = ddc::in_tags_v<typename FirstPatch::Grid1, Grid1DSeq>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange2,
                typename FirstPatch::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first(m_idx_ranges.template get<FirstPatch>());

        ddc::for_each(
                idx_range_par_first,
                [&](typename IdxRangeParFirstType::discrete_element_type const& idx_par) {
                    solve<eval_deriv>(functions_and_derivs, idx_par);
                });
    }


    /**
     * @brief Solve the matrix system MS = C to determine all the interface derivatives in the Grid1D
     * direction at a given index.
     * 
     * It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces. 
     *
     * @tparam IdxPar Type for the given index. It can be on the first or the second grid of the first patch.
     *
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @param[in] idx_par Index of the line in the Grid1D direction where we compute all the interface derivatives.
     */
    template <class IdxPar>
    void solve_cross_deriv(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        solve<eval_cross_deriv>(functions_and_derivs, idx_par);
    }

    /**
     * @brief Solve the matrix system MS = C to determine all the interface derivatives in the Grid1D
     * direction at all the last and the first indices of the first patch.
     * 
     * It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces. 
     *
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     */
    void solve_cross_deriv(MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        constexpr bool is_first_idx_range_perp_on_dim1
                = ddc::in_tags_v<typename FirstPatch::Grid1, Grid1DSeq>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange2,
                typename FirstPatch::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first(m_idx_ranges.template get<FirstPatch>());

        solve<eval_cross_deriv>(functions_and_derivs, idx_range_par_first.front());
        solve<eval_cross_deriv>(functions_and_derivs, idx_range_par_first.back());
    }



private:
    template <typename eval_type, class IdxPar>
    void solve(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        // Check the input types.
        static_assert(
                (std::is_same_v<IdxPar, typename FirstPatch::Idx1>)
                        || (std::is_same_v<IdxPar, typename FirstPatch::Idx2>),
                "The given index has to be a 1D index defined on the first patch. See "
                "InterfaceDerivativeMatrix<...>::first_patch.");

        // Update m_vector.
        set_vector<eval_type>(
                idx_par,
                functions_and_derivs,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});

        // Update m_interface_derivatives.
        m_solver->apply(m_vector, m_interface_derivatives);

        // Update derivatives.
        update_derivatives<eval_type>(
                functions_and_derivs,
                idx_par,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
    }

    /// @brief Set the matrix (I - M).
    template <std::size_t... I>
    void set_matrix(std::integer_sequence<std::size_t, I...>)
    {
        (set_line_matrix<I>(), ...);
    }

    /// @brief Fill in the line of the matrix corresponding to the given interface.
    template <std::size_t I>
    void set_line_matrix()
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        constexpr bool is_same_orientation = std::is_same_v<InterfaceI, EquivalentInterfaceI>;

        // Check that the Ith element of the derivative calculators collection is on InterfaceI.
        // using DerivativeCalculator
        //         = std::decay_t<std::tuple_element_t<I, SingleInterfaceDerivativesCalculatorTuple>>;
        // using InterfaceDerivativeCalculator = typename DerivativeCalculator::associated_interface;
        // using InterfaceDerivativeCalculatorSymmetrical = Interface<
        //         typename InterfaceDerivativeCalculator::Edge2,
        //         typename InterfaceDerivativeCalculator::Edge1,
        //         InterfaceDerivativeCalculator::orientations_agree>;
        // static_assert(
        //         (is_same_orientation && std::is_same_v<InterfaceI, InterfaceDerivativeCalculator>)
        //                 || (!is_same_orientation
        //                     && std::is_same_v<
        //                             InterfaceI,
        //                             InterfaceDerivativeCalculatorSymmetrical>),
        //         "The list of interface derivative calculators is not well sorted.");

        double const coeff_deriv_patch1
                = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                          .get_coeff_deriv_patch_1();
        double const coeff_deriv_patch2
                = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                          .get_coeff_deriv_patch_2();
        const double coeff_left = is_same_orientation ? coeff_deriv_patch1 : coeff_deriv_patch2;
        const double coeff_right = is_same_orientation ? coeff_deriv_patch2 : coeff_deriv_patch1;
        const double coefs[] = {-coeff_left, 1, -coeff_right};

        for (int dofs : {-1, 0, 1}) {
            if (0 <= I + dofs && I + dofs < n_inner_interfaces) {
                m_matrix->at(I, I + dofs) = coefs[dofs + 1];
            }
        }

        if constexpr (is_periodic && I == 0) {
            m_matrix->at(I, n_inner_interfaces - 1) = coefs[0];
        }
        if constexpr (is_periodic && I == n_inner_interfaces - 1) {
            m_matrix->at(I, 0) = coefs[2];
        }
    }


    /// @brief Set the vector C.
    template <typename eval_type, class GridPar1D, std::size_t... I>
    void set_vector(
            Idx<GridPar1D> const& slice_idx,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            std::integer_sequence<std::size_t, I...>)
    {
        static_assert(
                (std::is_same_v<GridPar1D, typename FirstPatch::Grid1>)
                || (std::is_same_v<GridPar1D, typename FirstPatch::Grid2>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim1>>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim2>>));

        IdxRange<GridPar1D> idx_range_1d_first_patch(m_idx_ranges.template get<FirstPatch>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (set_line_vector<I, eval_type>(slice_idx_value, functions_and_derivs), ...);
    }

    /// @brief Fill in the Ith line of the vector C corresponding.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_deriv>, bool> = true>
    void set_line_vector(
            int& slice_idx_1_value,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation of the ordered interface and the one in the derivative calculator.
        constexpr bool is_same_orientation = std::is_same_v<InterfaceI, EquivalentInterfaceI>;

        // If the orientations are not the same, we change the sign.
        const int sign = is_same_orientation - !is_same_orientation;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        constexpr Extremity extermity_1 = InterfaceI::Edge1::extremity;
        constexpr Extremity extermity_2 = InterfaceI::Edge2::extremity;

        using GridPerp1 = typename InterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename InterfaceI::Edge2::perpendicular_grid;

        using DimPerp1 = typename GridPerp1::continuous_dimension_type;
        using DimPerp2 = typename GridPerp2::continuous_dimension_type;

        using GridPar1 = typename InterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename InterfaceI::Edge2::parallel_grid;

        using DerivPerp1 = typename ddc::Deriv<DimPerp1>;
        using DerivPerp2 = typename ddc::Deriv<DimPerp2>;

        auto [idx_slice_1, idx_slice_2] = get_slice_indexes<InterfaceI>(slice_idx_1_value);

        DerivFieldOnPatch_host<Patch1> function_and_derivs_1
                = functions_and_derivs.template get<Patch1>();
        DerivFieldOnPatch_host<Patch2> function_and_derivs_2
                = functions_and_derivs.template get<Patch2>();


        // Use the function values to compute the first derivatives.
        DField<typename Patch1::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_1
                = function_and_derivs_1.get_values_field();
        DField<typename Patch2::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_2
                = function_and_derivs_2.get_values_field();

        double const lin_comb_funct
                = sign
                  * m_derivatives_calculators.template get<EquivalentInterfaceI>()
                            .get_function_coefficients(
                                    get_const_field(function_1[idx_slice_1]),
                                    get_const_field(function_2[idx_slice_2]));

        m_vector->get_values()[I] = lin_comb_funct;

        // TODO: Test this part.
        constexpr bool is_lower_bound_deriv_dependent
                = (LowerBound == ddc::BoundCond::HERMITE && I == 0);
        constexpr bool is_upper_bound_deriv_dependent
                = (UpperBound == ddc::BoundCond::HERMITE && I == n_inner_interfaces - 1);

        if constexpr (is_lower_bound_deriv_dependent || is_upper_bound_deriv_dependent) {
            IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch1>());
            IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch2>());

            IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1(
                    idx_range_perp_1.front(),
                    IdxStep<GridPerp1>(2),
                    idx_range_perp_1.extents());

            IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2(
                    idx_range_perp_2.front(),
                    IdxStep<GridPerp2>(2),
                    idx_range_perp_2.extents());

            Idx<GridPerp1> idx_deriv_1 = (extermity_1 == Extremity::BACK)
                                                 ? idx_range_slice_dperp_1.front()
                                                 : idx_range_slice_dperp_1.back();
            Idx<GridPerp2> idx_deriv_2 = (extermity_2 == Extremity::BACK)
                                                 ? idx_range_slice_dperp_2.front()
                                                 : idx_range_slice_dperp_2.back();

            Idx<DerivPerp1, GridPerp1> idx_slice_deriv_1(Idx<DerivPerp1>(1), idx_deriv_1);
            Idx<DerivPerp2, GridPerp2> idx_slice_deriv_2(Idx<DerivPerp2>(1), idx_deriv_2);

            const double deriv_1 = function_and_derivs_1[idx_slice_deriv_1](idx_slice_1);
            const double deriv_2 = function_and_derivs_2[idx_slice_deriv_2](idx_slice_2);

            const double coeff_deriv_1
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_1();
            const double coeff_deriv_2
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_2();

            if constexpr (is_lower_bound_deriv_dependent) {
                m_vector->get_values()[I]
                        += is_same_orientation ? coeff_deriv_1 * deriv_1 : -coeff_deriv_2 * deriv_2;
            } else {
                m_vector->get_values()[I]
                        += is_same_orientation ? coeff_deriv_2 * deriv_2 : -coeff_deriv_1 * deriv_1;
            }
        }

        //     if constexpr (LowerBound == ddc::BoundCond::HERMITE && I == 0) {
        //         {
        //             // auto [deriv_1, deriv_2]
        //             //         = get_interface_derivatives<I>(derivs_min, derivs_min, Extremity::BACK);
        //             // auto [idx_first_deriv_1, idx_first_deriv_2]
        //             //         = get_deriv_indexes<I>(deriv_1, deriv_2);

        //             IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch1>());
        //             IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch2>());

        //             IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1(
        //                     idx_range_perp_1.front(),
        //                     IdxStep<GridPerp1>(2),
        //                     idx_range_perp_1.extents());

        //             IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2(
        //                     idx_range_perp_1.front(),
        //                     IdxStep<GridPerp2>(2),
        //                     idx_range_perp_1.extents());

        //             Idx<GridPerp1> idx_deriv_1 = (extermity_1 == Extremity::BACK)
        //                                                  ? idx_range_slice_dperp_1.front()
        //                                                  : idx_range_slice_dperp_1.back();
        //             Idx<GridPerp2> idx_deriv_2 = (extermity_2 == Extremity::BACK)
        //                                                  ? idx_range_slice_dperp_2.front()
        //                                                  : idx_range_slice_dperp_2.back();

        //             Idx<DerivPerp1, GridPerp1> idx_slice_deriv_1(Idx<DerivPerp1>(1), idx_deriv_1);
        //             Idx<DerivPerp2, GridPerp2> idx_slice_deriv_2(Idx<DerivPerp2>(1), idx_deriv_2);

        //             // DField<IdxRange<GridPar1>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_1
        //             //         = function_and_derivs_1[idx_slice_deriv_1];
        //             // DField<IdxRange<GridPar2>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_2
        //             //         = function_and_derivs_2[idx_slice_deriv_2];

        //             const double deriv_1 = function_and_derivs_1[idx_slice_deriv_1](idx_slice_1);
        //             const double deriv_2 = function_and_derivs_2[idx_slice_deriv_2](idx_slice_2);

        //             double coeff_deriv_1
        //                     = m_derivatives_calculators.template get<EquivalentInterfaceI>()
        //                               .get_coeff_deriv_patch_1();
        //             double coeff_deriv_2
        //                     = m_derivatives_calculators.template get<EquivalentInterfaceI>()
        //                               .get_coeff_deriv_patch_2();

        //             m_vector->get_values()[I]
        //                     += is_same_orientation ? coeff_deriv_1 * deriv_1 : -coeff_deriv_2 * deriv_2;
        //         }
        //     }
        //     if constexpr (UpperBound == ddc::BoundCond::HERMITE && I == n_inner_interfaces - 1) {
        //         // auto [deriv_1, deriv_2]
        //         //         = get_interface_derivatives<I>(derivs_min, derivs_min, Extremity::BACK);
        //         // auto [idx_first_deriv_1, idx_first_deriv_2] = get_deriv_indexes<I>(deriv_1, deriv_2);

        //         // double coeff_deriv_1 = m_derivatives_calculators.template get<EquivalentInterfaceI>()
        //         //                                .get_coeff_deriv_patch_1();
        //         // double coeff_deriv_2 = m_derivatives_calculators.template get<EquivalentInterfaceI>()
        //         //                                .get_coeff_deriv_patch_2();

        //         // m_vector->get_values()[I]
        //         //         += is_same_orientation
        //         //                    ? coeff_deriv_2 * deriv_2(idx_first_deriv_2, idx_slice_2)
        //         //                    : -coeff_deriv_1 * deriv_1(idx_first_deriv_1, idx_slice_1);
        //         IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch1>());
        //         IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch2>());

        //         IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1(
        //                 idx_range_perp_1.front(),
        //                 IdxStep<GridPerp1>(2),
        //                 idx_range_perp_1.extents());

        //         IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2(
        //                 idx_range_perp_1.front(),
        //                 IdxStep<GridPerp2>(2),
        //                 idx_range_perp_1.extents());

        //         Idx<GridPerp1> idx_deriv_1 = (extermity_1 == Extremity::BACK)
        //                                              ? idx_range_slice_dperp_1.front()
        //                                              : idx_range_slice_dperp_1.back();
        //         Idx<GridPerp2> idx_deriv_2 = (extermity_2 == Extremity::BACK)
        //                                              ? idx_range_slice_dperp_2.front()
        //                                              : idx_range_slice_dperp_2.back();

        //         Idx<DerivPerp1, GridPerp1> idx_slice_deriv_1(Idx<DerivPerp1>(1), idx_deriv_1);
        //         Idx<DerivPerp2, GridPerp2> idx_slice_deriv_2(Idx<DerivPerp2>(1), idx_deriv_2);

        //         // DField<IdxRange<GridPar1>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_1
        //         //         = function_and_derivs_1[idx_slice_deriv_1];
        //         // DField<IdxRange<GridPar2>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_2
        //         //         = function_and_derivs_2[idx_slice_deriv_2];

        //         const double deriv_1 = function_and_derivs_1[idx_slice_deriv_1](idx_slice_1);
        //         const double deriv_2 = function_and_derivs_2[idx_slice_deriv_2](idx_slice_2);

        //         double coeff_deriv_1 = m_derivatives_calculators.template get<EquivalentInterfaceI>()
        //                                        .get_coeff_deriv_patch_1();
        //         double coeff_deriv_2 = m_derivatives_calculators.template get<EquivalentInterfaceI>()
        //                                        .get_coeff_deriv_patch_2();

        //         m_vector->get_values()[I]
        //                 += is_same_orientation ? coeff_deriv_2 * deriv_2 : -coeff_deriv_1 * deriv_1;
        //     }
    }

    /// @brief Fill in the Ith line of the vector C corresponding.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_cross_deriv>, bool> = true>
    void set_line_vector(
            int& slice_idx_1_value,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation of the ordered interface and the one in the derivative calculator.
        constexpr bool is_same_orientation = std::is_same_v<InterfaceI, EquivalentInterfaceI>;

        // If the orientations are not the same, we change the sign.
        const int sign = is_same_orientation - !is_same_orientation;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        auto [idx_slice_1, idx_slice_2] = get_slice_indexes<InterfaceI>(slice_idx_1_value);



        DerivFieldOnPatch_host<Patch1> function_and_derivs_1
                = functions_and_derivs.template get<Patch1>();
        DerivFieldOnPatch_host<Patch2> function_and_derivs_2
                = functions_and_derivs.template get<Patch2>();

        // Use the first derivatives to compute the cross-derivatives.
        using GridPerp1 = typename InterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename InterfaceI::Edge2::perpendicular_grid;

        using GridPar1 = typename InterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename InterfaceI::Edge2::parallel_grid;

        using DimPar1 = typename GridPar1::continuous_dimension_type;
        using DimPar2 = typename GridPar2::continuous_dimension_type;

        using DerivPar1 = typename ddc::Deriv<DimPar1>;
        using DerivPar2 = typename ddc::Deriv<DimPar2>;

        IdxRange<GridPar1> idx_range_par_1(m_idx_ranges.template get<Patch1>());
        IdxRange<GridPar2> idx_range_par_2(m_idx_ranges.template get<Patch2>());

        // The slice indices has to be a point at a corner, i.e. index range boundaries.
        assert((idx_slice_1 == idx_range_par_1.front()) || (idx_slice_1 == idx_range_par_1.back()));
        assert((idx_slice_2 == idx_range_par_2.front()) || (idx_slice_2 == idx_range_par_2.back()));


        IdxRangeSlice<GridPar1> idx_range_slice_dpar_1(
                idx_range_par_1.front(),
                IdxStep<GridPar1>(2),
                idx_range_par_1.extents());

        IdxRangeSlice<GridPar2> idx_range_slice_dpar_2(
                idx_range_par_2.front(),
                IdxStep<GridPar2>(2),
                idx_range_par_2.extents());

        Idx<GridPar1> idx_deriv_par_1 = (idx_slice_1 == idx_range_par_1.front())
                                                ? idx_range_slice_dpar_1.front()
                                                : idx_range_slice_dpar_1.back();
        Idx<GridPar2> idx_deriv_par_2 = (idx_slice_2 == idx_range_par_2.front())
                                                ? idx_range_slice_dpar_2.front()
                                                : idx_range_slice_dpar_2.back();

        Idx<DerivPar1, GridPar1> idx_slice_deriv_1(Idx<DerivPar1>(1), idx_deriv_par_1);
        Idx<DerivPar2, GridPar2> idx_slice_deriv_2(Idx<DerivPar2>(1), idx_deriv_par_2);

        DField<IdxRange<GridPerp1>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_1
                = function_and_derivs_1[idx_slice_deriv_1];
        DField<IdxRange<GridPerp2>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_2
                = function_and_derivs_2[idx_slice_deriv_2];

        double const lin_comb_funct
                = sign
                  * m_derivatives_calculators.template get<EquivalentInterfaceI>()
                            .get_function_coefficients(
                                    get_const_field(derivs_1),
                                    get_const_field(derivs_2));

        m_vector->get_values()[I] = lin_comb_funct;

        // TODO: Test this part.
        constexpr bool is_lower_bound_deriv_dependent
                = (LowerBound == ddc::BoundCond::HERMITE && I == 0);
        constexpr bool is_upper_bound_deriv_dependent
                = (UpperBound == ddc::BoundCond::HERMITE && I == n_inner_interfaces - 1);

        if constexpr (is_lower_bound_deriv_dependent || is_upper_bound_deriv_dependent) {
            constexpr Extremity extermity_1 = InterfaceI::Edge1::extremity;
            constexpr Extremity extermity_2 = InterfaceI::Edge2::extremity;

            using Grid1_1 = typename Patch1::Grid1;
            using Grid2_1 = typename Patch1::Grid2;
            using Grid1_2 = typename Patch2::Grid1;
            using Grid2_2 = typename Patch2::Grid2;

            using Deriv1_1 = ddc::Deriv<typename Patch1::Dim1>;
            using Deriv2_1 = ddc::Deriv<typename Patch1::Dim2>;
            using Deriv1_2 = ddc::Deriv<typename Patch2::Dim1>;
            using Deriv2_2 = ddc::Deriv<typename Patch2::Dim2>;

            constexpr bool is_grid_par_1_on_dim1 = std::is_same_v<GridPar1, Grid1_1>;
            constexpr bool is_grid_par_2_on_dim1 = std::is_same_v<GridPar2, Grid1_2>;

            IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch1>());
            IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch2>());

            IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1(
                    idx_range_perp_1.front(),
                    IdxStep<GridPerp1>(2),
                    idx_range_perp_1.extents());

            IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2(
                    idx_range_perp_2.front(),
                    IdxStep<GridPerp2>(2),
                    idx_range_perp_2.extents());

            Idx<GridPerp1> idx_deriv_perp_1 = (extermity_1 == Extremity::BACK)
                                                      ? idx_range_slice_dperp_1.front()
                                                      : idx_range_slice_dperp_1.back();
            Idx<GridPerp2> idx_deriv_perp_2 = (extermity_2 == Extremity::BACK)
                                                      ? idx_range_slice_dperp_2.front()
                                                      : idx_range_slice_dperp_2.back();

            Idx<Grid1_1, Grid2_1> idx_d1d2_1;
            if constexpr (is_grid_par_1_on_dim1) {
                idx_d1d2_1 = Idx<Grid1_1, Grid2_1>(idx_deriv_par_1, idx_deriv_perp_1);
            } else {
                idx_d1d2_1 = Idx<Grid1_1, Grid2_1>(idx_deriv_perp_1, idx_deriv_par_1);
            }

            Idx<Grid1_2, Grid2_2> idx_d1d2_2;
            if constexpr (is_grid_par_2_on_dim1) {
                idx_d1d2_2 = Idx<Grid1_2, Grid2_2>(idx_deriv_par_2, idx_deriv_perp_2);
            } else {
                idx_d1d2_2 = Idx<Grid1_2, Grid2_2>(idx_deriv_perp_2, idx_deriv_par_2);
            }

            Idx<Grid1_1> idx_d1_1(idx_d1d2_1);
            Idx<Grid2_1> idx_d2_1(idx_d1d2_1);
            Idx<Grid1_2> idx_d1_2(idx_d1d2_2);
            Idx<Grid2_2> idx_d2_2(idx_d1d2_2);

            Idx<Deriv1_1, Grid1_1, Deriv2_1, Grid2_1>
                    idx_cross_deriv1(Idx<Deriv1_1>(1), idx_d1_1, Idx<Deriv2_1>(1), idx_d2_1);
            Idx<Deriv1_2, Grid1_2, Deriv2_2, Grid2_2>
                    idx_cross_deriv2(Idx<Deriv1_2>(1), idx_d1_2, Idx<Deriv2_2>(1), idx_d2_2);

            const double cross_deriv_1 = function_and_derivs_1(idx_cross_deriv1);
            const double cross_deriv_2 = function_and_derivs_2(idx_cross_deriv2);

            const double coeff_deriv_1
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_1();
            const double coeff_deriv_2
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_2();

            if constexpr (is_lower_bound_deriv_dependent) {
                m_vector->get_values()[I] += is_same_orientation ? coeff_deriv_1 * cross_deriv_1
                                                                 : -coeff_deriv_2 * cross_deriv_2;
            } else {
                m_vector->get_values()[I] += is_same_orientation ? coeff_deriv_2 * cross_deriv_2
                                                                 : -coeff_deriv_1 * cross_deriv_1;
            }
        }
    }


    /// @brief Once the interface derivatives computed, fill in the correct derivative fields.
    template <typename eval_type, class GridPar1D, std::size_t... I>
    void update_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            Idx<GridPar1D> const& slice_idx,
            std::integer_sequence<std::size_t, I...>)
    {
        static_assert(
                (std::is_same_v<GridPar1D, typename FirstPatch::Grid1>)
                || (std::is_same_v<GridPar1D, typename FirstPatch::Grid2>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim1>>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim2>>));

        IdxRange<GridPar1D> idx_range_1d_first_patch(m_idx_ranges.template get<FirstPatch>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (update_derivatives_at_interface<I, eval_type>(functions_and_derivs, slice_idx_value), ...);
    }


    /// @brief Associate the Ith derivative values to the correct first derivative.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_deriv>, bool> = true>
    void update_derivatives_at_interface(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& slice_idx_1_value)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        constexpr Extremity extermity_1 = InterfaceI::Edge1::extremity;
        constexpr Extremity extermity_2 = InterfaceI::Edge2::extremity;

        using GridPerp1 = typename InterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename InterfaceI::Edge2::perpendicular_grid;

        using GridPar1 = typename InterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename InterfaceI::Edge2::parallel_grid;

        using DimPerp1 = typename GridPerp1::continuous_dimension_type;
        using DimPerp2 = typename GridPerp2::continuous_dimension_type;

        using DerivPerp1 = typename ddc::Deriv<DimPerp1>;
        using DerivPerp2 = typename ddc::Deriv<DimPerp2>;

        // Get the fields of the left and right patch of the interface.
        DerivFieldOnPatch_host<Patch1> function_and_derivs_1
                = functions_and_derivs.template get<Patch1>();
        DerivFieldOnPatch_host<Patch2> function_and_derivs_2
                = functions_and_derivs.template get<Patch2>();

        // Get the correct index ranges and indices for the slices.
        IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch1>());
        IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch2>());

        IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1(
                idx_range_perp_1.front(),
                IdxStep<GridPerp1>(2),
                idx_range_perp_1.extents());

        IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2(
                idx_range_perp_2.front(),
                IdxStep<GridPerp2>(2),
                idx_range_perp_2.extents());

        Idx<GridPerp1> idx_deriv_1 = (extermity_1 == Extremity::FRONT)
                                             ? idx_range_slice_dperp_1.front()
                                             : idx_range_slice_dperp_1.back();
        Idx<GridPerp2> idx_deriv_2 = (extermity_2 == Extremity::FRONT)
                                             ? idx_range_slice_dperp_2.front()
                                             : idx_range_slice_dperp_2.back();

        Idx<DerivPerp1, GridPerp1> idx_slice_deriv_1(Idx<DerivPerp1>(1), idx_deriv_1);
        Idx<DerivPerp2, GridPerp2> idx_slice_deriv_2(Idx<DerivPerp2>(1), idx_deriv_2);

        auto [slice_idx_1, slice_idx_2] = get_slice_indexes<InterfaceI>(slice_idx_1_value);

        // Select the correct data and update.
        DField<IdxRange<GridPar1>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_1
                = function_and_derivs_1[idx_slice_deriv_1];
        DField<IdxRange<GridPar2>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_2
                = function_and_derivs_2[idx_slice_deriv_2];

        deriv_1(slice_idx_1) = m_interface_derivatives->get_values()[I];
        deriv_2(slice_idx_2) = m_interface_derivatives->get_values()[I];
    }

    /// @brief Associate the Ith derivative values to the correct cross-derivative.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_cross_deriv>, bool> = true>
    void update_derivatives_at_interface(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& slice_idx_1_value)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        constexpr Extremity extermity_1 = InterfaceI::Edge1::extremity;
        constexpr Extremity extermity_2 = InterfaceI::Edge2::extremity;

        using GridPar1 = typename InterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename InterfaceI::Edge2::parallel_grid;

        using Grid1_1 = typename Patch1::Grid1;
        using Grid2_1 = typename Patch1::Grid2;
        using Grid1_2 = typename Patch2::Grid1;
        using Grid2_2 = typename Patch2::Grid2;

        using Deriv1_1 = ddc::Deriv<typename Patch1::Dim1>;
        using Deriv2_1 = ddc::Deriv<typename Patch1::Dim2>;
        using Deriv1_2 = ddc::Deriv<typename Patch2::Dim1>;
        using Deriv2_2 = ddc::Deriv<typename Patch2::Dim2>;

        constexpr bool is_grid_par_1_on_dim1 = std::is_same_v<GridPar1, Grid1_1>;
        constexpr bool is_grid_par_2_on_dim1 = std::is_same_v<GridPar2, Grid1_2>;

        // Get the fields of the left and right patch of the interface.
        DerivFieldOnPatch_host<Patch1> function_and_derivs_1
                = functions_and_derivs.template get<Patch1>();
        DerivFieldOnPatch_host<Patch2> function_and_derivs_2
                = functions_and_derivs.template get<Patch2>();

        // Get the correct index ranges and indices for the slices.
        IdxRange<Grid1_1> idx_range_grid1_1(m_idx_ranges.template get<Patch1>());
        IdxRange<Grid2_1> idx_range_grid2_1(m_idx_ranges.template get<Patch1>());
        IdxRange<Grid1_2> idx_range_grid1_2(m_idx_ranges.template get<Patch2>());
        IdxRange<Grid2_2> idx_range_grid2_2(m_idx_ranges.template get<Patch2>());

        IdxRangeSlice<Grid1_1> idx_range_slice_d1_1(
                idx_range_grid1_1.front(),
                IdxStep<Grid1_1>(2),
                idx_range_grid1_1.extents());
        IdxRangeSlice<Grid2_1> idx_range_slice_d2_1(
                idx_range_grid2_1.front(),
                IdxStep<Grid2_1>(2),
                idx_range_grid2_1.extents());

        IdxRangeSlice<Grid1_2> idx_range_slice_d1_2(
                idx_range_grid1_2.front(),
                IdxStep<Grid1_2>(2),
                idx_range_grid1_2.extents());
        IdxRangeSlice<Grid2_2> idx_range_slice_d2_2(
                idx_range_grid2_2.front(),
                IdxStep<Grid2_2>(2),
                idx_range_grid2_2.extents());

        IdxRange<GridPar1> idx_range_par_1(m_idx_ranges.template get<Patch1>());
        IdxRange<GridPar2> idx_range_par_2(m_idx_ranges.template get<Patch2>());

        auto [slice_idx_1, slice_idx_2] = get_slice_indexes<InterfaceI>(slice_idx_1_value);

        // The slice indices has to be a point at a corner, i.e. index range boundaries.
        assert((slice_idx_1 == idx_range_par_1.front()) || (slice_idx_1 == idx_range_par_1.back()));
        assert((slice_idx_2 == idx_range_par_2.front()) || (slice_idx_2 == idx_range_par_2.back()));

        const bool is_idx_par_min_1 = (slice_idx_1 == idx_range_par_1.front());
        const bool is_idx_par_min_2 = (slice_idx_2 == idx_range_par_2.front());

        // --- Update the cross-derivative on Patch1
        Idx<Grid1_1> idx_d1_1;
        Idx<Grid2_1> idx_d2_1;
        if (is_idx_par_min_1 && (extermity_1 == Extremity::FRONT)) {
            // corner min min
            idx_d1_1 = idx_range_slice_d1_1.front();
            idx_d2_1 = idx_range_slice_d2_1.front();
        } else if (!is_idx_par_min_1 && (extermity_1 == Extremity::BACK)) {
            // corner max max
            idx_d1_1 = idx_range_slice_d1_1.back();
            idx_d2_1 = idx_range_slice_d2_1.back();
        } else if (
                (is_grid_par_1_on_dim1 && (extermity_1 == Extremity::BACK))
                || (!is_grid_par_1_on_dim1 && (extermity_1 == Extremity::FRONT))) {
            // corner min max
            idx_d1_1 = idx_range_slice_d1_1.front();
            idx_d2_1 = idx_range_slice_d2_1.back();
        } else {
            // corner max min
            idx_d1_1 = idx_range_slice_d1_1.back();
            idx_d2_1 = idx_range_slice_d2_1.front();
        }
        Idx<Deriv1_1, Grid1_1, Deriv2_1, Grid2_1>
                idx_cross_deriv1(Idx<Deriv1_1>(1), idx_d1_1, Idx<Deriv2_1>(1), idx_d2_1);
        function_and_derivs_1(idx_cross_deriv1) = m_interface_derivatives->get_values()[I];

        // --- Update the cross-derivative on Patch2
        Idx<Grid1_2> idx_d1_2;
        Idx<Grid2_2> idx_d2_2;
        if (is_idx_par_min_2 && (extermity_2 == Extremity::FRONT)) {
            // corner min min
            idx_d1_2 = idx_range_slice_d1_2.front();
            idx_d2_2 = idx_range_slice_d2_2.front();
        } else if (!is_idx_par_min_2 && (extermity_2 == Extremity::BACK)) {
            // corner max max
            idx_d1_2 = idx_range_slice_d1_2.back();
            idx_d2_2 = idx_range_slice_d2_2.back();
        } else if (
                (is_grid_par_2_on_dim1 && (extermity_2 == Extremity::BACK))
                || (!is_grid_par_2_on_dim1 && (extermity_2 == Extremity::FRONT))) {
            // corner min max
            idx_d1_2 = idx_range_slice_d1_2.front();
            idx_d2_2 = idx_range_slice_d2_2.back();
        } else {
            // corner max min
            idx_d1_2 = idx_range_slice_d1_2.back();
            idx_d2_2 = idx_range_slice_d2_2.front();
        }
        Idx<Deriv1_2, Grid1_2, Deriv2_2, Grid2_2>
                idx_cross_deriv2(Idx<Deriv1_2>(1), idx_d1_2, Idx<Deriv2_2>(1), idx_d2_2);
        function_and_derivs_2(idx_cross_deriv2) = m_interface_derivatives->get_values()[I];
    }


    /// @brief Get the first derivative field on Patch1 and Patch2 for the Ith interface.
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


    /// @brief @brief Get the Patch1 and Patch2 derivative indices on the grid parallel to the interface.
    template <std::size_t I, class DerivFieldOnPatch1, class DerivFieldOnPatch2>
    auto get_deriv_indexes(DerivFieldOnPatch1 const& deriv_1, DerivFieldOnPatch2 const& deriv_2)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        using PerpGrid1 = typename InterfaceI::Edge1::perpendicular_grid;
        using PerpGrid2 = typename InterfaceI::Edge2::perpendicular_grid;

        using PerpDim1 = typename PerpGrid1::continuous_dimension_type;
        using PerpDim2 = typename PerpGrid2::continuous_dimension_type;

        // Determine the parallel indexes.
        IdxRange<ddc::Deriv<PerpDim1>> idx_range_first_deriv_1(get_idx_range(deriv_1));
        IdxRange<ddc::Deriv<PerpDim2>> idx_range_first_deriv_2(get_idx_range(deriv_2));

        Idx<ddc::Deriv<PerpDim1>> idx_first_deriv_1(idx_range_first_deriv_1.front());
        Idx<ddc::Deriv<PerpDim2>> idx_first_deriv_2(idx_range_first_deriv_2.front());

        return std::make_tuple(idx_first_deriv_1, idx_first_deriv_2);
    }

    template <typename InterfaceI>
    auto get_slice_indexes(int& slice_idx_1_value)
    {
        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        using ParallGrid1 = typename InterfaceI::Edge1::parallel_grid;
        using ParallGrid2 = typename InterfaceI::Edge2::parallel_grid;

        IdxRange<ParallGrid1> const idx_range_parell_1(m_idx_ranges.template get<Patch1>());
        IdxRange<ParallGrid2> const idx_range_parell_2(m_idx_ranges.template get<Patch2>());

        EdgeTransformation<InterfaceI> index_converter(idx_range_parell_1, idx_range_parell_2);

        Idx<ParallGrid1> slice_idx_1(slice_idx_1_value);
        Idx<ParallGrid2> slice_idx_2(index_converter(slice_idx_1));
        slice_idx_1_value = (slice_idx_2 - idx_range_parell_2.front()).value();
        return std::make_tuple(slice_idx_1, slice_idx_2);
    }
};
