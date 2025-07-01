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
        template <typename P>
        typename ValuesOnPatch, // Fix type ? 
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
            ((LowerBound == ddc::BoundCond::PERIODIC) && (UpperBound == ddc::BoundCond::PERIODIC))
                    || ((LowerBound != ddc::BoundCond::PERIODIC)
                        && (UpperBound != ddc::BoundCond::PERIODIC)),
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
            "The number of 1D grids and inner interfaces should fit.");

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

    // Grid perpendicular to the interface for a given patch.
    template <typename Patch>
    using GridPerpOnPatch = typename GetDerivsOnPatch<Patch>::GridPerp;

    // Grid parallel to the interface for a given patch.
    template <typename Patch>
    using GridParOnPatch = typename GetDerivsOnPatch<Patch>::GridPar;

    // Field for the first derivatives along the perpendicular direction of the interface.
    template <typename Patch>
    using DerivsPerpOnPatch = typename GetDerivsOnPatch<Patch>::type_perp;

    // Const field the first derivatives along the parallel direction of the interface.
    template <typename Patch>
    using ConstDerivsParOnPatch = typename GetDerivsOnPatch<Patch>::type_par;

    // Field for the cross-derivatives.
    template <typename Patch>
    using CrossDerivsPerpOnPatch = DField<
            IdxRange<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>>,
            typename ExecSpace::memory_space>;


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

    // TODO: WARNING, BE ABLE TO DEAL WITH DIFFERENT ORIENTATION INTERFACES.
    using SingleInterfaceDerivativesCalculatorTuple
            = get_tuple_deriv_calculator_t<inner_interface_collection>;
    SingleInterfaceDerivativesCalculatorTuple const& m_derivatives_calculators;

public:
    ~InterfaceDerivativeMatrix() {}

    /**
     * @brief Instantiate InterfaceDerivativeMatrix. 
     * 
     * It creates a dense Gingko matrix, a vector for the rhs and a vector for the solution. 
     * The matrix is initialised with the derivatives calculator collection given in input. 
     * 
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param derivatives_calculators Tuple of SingleInterfaceDerivativesCalculator containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     * @warning The interface derivative calculators have to be ordered in the tuple. 
     *      The order has to be the same as the interfaces in the Grid1D direction. 
     */
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


    /**
     * @brief Solve the matrix system MS = C to determine all the interface derivatives in the Grid1D 
     * direction at a given index. 
     * 
     * @tparam IdxPar Type for the given index. It can be on the first or the second grid of the first patch. 
     * @tparam FieldOnPatch Templated on the patch, it is the type for constant field of function values 
     *              or first derivatives on the patch list. 
     * @tparam DerivsTypeCollection MultipatchField or std::tuple of first derivatives or cross-derivatives 
     *              we want to compute.
     * 
     * @param derivs_min[out] Collection first derivatives or cross-derivatives on the west/south side we want to compute.
     * @param derivs_max[out] Collection first derivatives or cross-derivatives on the east/north side we want to compute.
     * @param idx_par[in] Index of the line in the Grid1D direction where we compute all the interface derivatives. 
     * @param function_values[in] Constant field collection of function values or first derivatives on the patch list. 
     */
    template <class IdxPar, template <typename P> class FieldOnPatch, class DerivsTypeCollection>
    void solve( // Should be useful for non-conforming case (later).
            DerivsTypeCollection const& derivs_min,
            DerivsTypeCollection const& derivs_max,
            IdxPar const& idx_par,
            MultipatchField<FieldOnPatch, Patches...> const& function_values)
    {
        // Check the input types.
        static_assert(
                (std::is_same_v<IdxPar, typename FirstPatch::Idx1>)
                        || (std::is_same_v<IdxPar, typename FirstPatch::Idx2>),
                "The given index has to be a 1D index defined on the first patch. See "
                "InterfaceDerivativeMatrix<...>::first_patch.");
        static_assert(
                (std::is_same_v<FieldOnPatch<FirstPatch>, ValuesOnPatch<FirstPatch>>)
                        || (std::is_same_v<
                                FieldOnPatch<FirstPatch>,
                                ConstDerivsParOnPatch<FirstPatch>>),
                "The `function_values` parameter has to be a MultipatchField collection of "
                "function values (see ValuesOnPatch) or first derivative fields.");
        // RMK: It seems to not recognise MultipatchField<DerivsPerpOnPatch, Patches...> being equal to
        // MultipatchField<Deriv[1/2]_OnPatch_2D_host, Patches...>
        static_assert(
                (std::is_same_v<
                        DerivsTypeCollection,
                        MultipatchField<DerivsPerpOnPatch, Patches...>>)
                        || (std::is_same_v<
                                DerivsTypeCollection,
                                MultipatchField<Deriv1_OnPatch_2D_host, Patches...>>)
                        || (std::is_same_v<
                                DerivsTypeCollection,
                                MultipatchField<Deriv2_OnPatch_2D_host, Patches...>>)
                        || (std::is_same_v<
                                DerivsTypeCollection,
                                MultipatchField<Deriv1_OnPatch_2D, Patches...>>)
                        || (std::is_same_v<
                                DerivsTypeCollection,
                                MultipatchField<Deriv2_OnPatch_2D, Patches...>>)
                        || (std::is_same_v<
                                DerivsTypeCollection,
                                std::tuple<CrossDerivsPerpOnPatch<Patches>...>>),
                "The `derivs_min/max` parameters are either first derivative field MultipatchField "
                "collection, or cross-derivative field tuple collection.");

        // Update m_vector.
        set_vector(
                idx_par,
                function_values,
                derivs_min,
                derivs_max,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});

        // Update m_interface_derivatives.
        m_solver->apply(m_vector, m_interface_derivatives);

        // Update derivatives.
        update_derivatives(
                derivs_min,
                derivs_max,
                function_values,
                idx_par,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
    }

    /**
     * @brief Solve the matrix system MS = C to determine all the interface derivatives in the Grid1D 
     * direction at all the indices of the first patch. 
     * 
     * @tparam FieldOnPatch Templated on the patch, it is the type for constant field of function values 
     *              or first derivatives on the patch list. 
     * @tparam DerivsTypeCollection MultipatchField or std::tuple of first derivatives or cross-derivatives 
     *              we want to compute.
     * 
     * @param derivs_min[out] Collection first derivatives or cross-derivatives on the west/south side we want to compute.
     * @param derivs_max[out] Collection first derivatives or cross-derivatives on the east/north side we want to compute.
     * @param function_values[in] Constant field collection of function values or first derivatives on the patch list. 
     */
    template <template <typename P> class FieldOnPatch, class DerivsTypeCollection>
    void solve(
            DerivsTypeCollection const& derivs_min,
            DerivsTypeCollection const& derivs_max,
            MultipatchField<FieldOnPatch, Patches...> const& function_values)
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
                    (*this).solve(derivs_min, derivs_max, idx_par, function_values);
                });
    }


private:
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

        constexpr bool is_same_orientation = std::is_same_v<
                InterfaceI,
                find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>>;

        // Check that the Ith element of the derivative calculators collection is on InterfaceI.
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
                            && std::is_same_v<
                                    InterfaceI,
                                    InterfaceDerivativeCalculatorSymmetrical>),
                "The list of interface derivative calculators is not well sorted.");

        const double coeff_left
                = is_same_orientation
                          ? std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1()
                          : std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();
        const double coeff_right
                = is_same_orientation
                          ? std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2()
                          : std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
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


    /// @brief Fill in the Ith line of the vector C corresponding.
    template <std::size_t I, template <typename P> class FieldOnPatch, class DerivCollectionType>
    void set_line_vector(
            int& slice_idx_1_value,
            MultipatchField<FieldOnPatch, Patches...> const& functions,
            DerivCollectionType const& derivs_min,
            DerivCollectionType const& derivs_max)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        // The orientation of the ordered interface and the one in the derivative calculator.
        constexpr bool is_same_orientation = std::is_same_v<
                InterfaceI,
                find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>>;

        // If the orientations are not the same, we change the sign.
        const int sign = is_same_orientation - !is_same_orientation;

        using Patch1 = typename InterfaceI::Edge1::associated_patch;
        using Patch2 = typename InterfaceI::Edge2::associated_patch;

        auto [idx_slice_1, idx_slice_2]
                = get_slice_indexes<I, Patch1, Patch2>(slice_idx_1_value, functions);

        FieldOnPatch<Patch1> function_1 = functions.template get<Patch1>();
        FieldOnPatch<Patch2> function_2 = functions.template get<Patch2>();

        double lin_comb_funct = sign
                                * std::get<I>(m_derivatives_calculators)
                                          .get_function_coefficients(
                                                  function_1[idx_slice_1],
                                                  function_2[idx_slice_2]);

        m_vector->get_values()[I] = lin_comb_funct;

        // TODO: Test this part.
        if constexpr (LowerBound == ddc::BoundCond::HERMITE && I == 0) {
            {
                auto [deriv_1, deriv_2]
                        = get_interface_derivatives<I>(derivs_min, derivs_min, Extremity::BACK);
                auto [idx_first_deriv_1, idx_first_deriv_2]
                        = get_deriv_indexes<I>(deriv_1, deriv_2);

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
            auto [idx_first_deriv_1, idx_first_deriv_2] = get_deriv_indexes<I>(deriv_1, deriv_2);

            double coeff_deriv_1 = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_1();
            double coeff_deriv_2 = std::get<I>(m_derivatives_calculators).get_coeff_deriv_patch_2();

            m_vector->get_values()[I]
                    += is_same_orientation
                               ? coeff_deriv_2 * deriv_2(idx_first_deriv_2, idx_slice_2)
                               : -coeff_deriv_1 * deriv_1(idx_first_deriv_1, idx_slice_1);
        }
    }

    /// @brief Once the interface derivatives computed, fill in the correct derivative fields.
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

    /// @brief Associate the Ith derivative values to the correct derivative field.
    template <std::size_t I, template <typename P> class FieldOnPatch, class DerivCollectionType>
    void update_derivatives_at_interface(
            DerivCollectionType const& derivs_min,
            DerivCollectionType const& derivs_max,
            MultipatchField<FieldOnPatch, Patches...> const& function_values,
            int& slice_idx_1_value)
    {
        auto [deriv_1, deriv_2]
                = get_interface_derivatives<I>(derivs_min, derivs_max, Extremity::FRONT);

        auto [slice_idx_1, slice_idx_2] = get_slice_indexes<I>(slice_idx_1_value, function_values);
        auto [idx_first_deriv_1, idx_first_deriv_2] = get_deriv_indexes<I>(deriv_1, deriv_2);

        // Update derivatives
        deriv_1(idx_first_deriv_1, slice_idx_1) = m_interface_derivatives->get_values()[I];
        deriv_2(idx_first_deriv_2, slice_idx_2) = m_interface_derivatives->get_values()[I];
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


    /// @brief Get the cross-derivative field on Patch1 and Patch2 for the Ith interface.
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

    /// @brief Get the slice indices on Patch1 and Patch2 corresponding to the given integer on Patch1.
    template <
            std::size_t I,
            class Patch1 = typename ddc::type_seq_element_t<I, inner_interface_collection>::Edge1::
                    associated_patch,
            class Patch2 = typename ddc::type_seq_element_t<I, inner_interface_collection>::Edge2::
                    associated_patch,
            template <typename P>
            class FieldOnPatch>
    auto get_slice_indexes(
            int& slice_idx_1_value,
            MultipatchField<FieldOnPatch, Patches...> const& function_values)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;

        using PerpGrid1 = typename InterfaceI::Edge1::perpendicular_grid;
        using PerpGrid2 = typename InterfaceI::Edge2::perpendicular_grid;

        using ParallGrid1 = typename InterfaceI::Edge1::parallel_grid;
        using ParallGrid2 = typename InterfaceI::Edge2::parallel_grid;

        // Define all the different index ranges.
        IdxRange<PerpGrid1> const idx_range_perp_1(m_idx_ranges.template get<Patch1>());
        IdxRange<PerpGrid2> const idx_range_perp_2(m_idx_ranges.template get<Patch2>());

        IdxRange<ParallGrid1> const idx_range_parell_1(m_idx_ranges.template get<Patch1>());
        IdxRange<ParallGrid2> const idx_range_parell_2(m_idx_ranges.template get<Patch2>());

        FieldOnPatch<Patch1> function_1 = function_values.template get<Patch1>();
        FieldOnPatch<Patch2> function_2 = function_values.template get<Patch2>();

        // Index range on Grid(s) or ddc::Deriv(s).
        auto const idx_range_fct_parall_1
                = ddc::remove_dims_of(get_idx_range(function_1), idx_range_perp_1);
        auto const idx_range_fct_parall_2
                = ddc::remove_dims_of(get_idx_range(function_2), idx_range_perp_2);

        // Get slice indexes.
        EdgeTransformation<InterfaceI> index_converter(idx_range_parell_1, idx_range_parell_2);

        using OIdx1 = typename decltype(idx_range_fct_parall_1)::discrete_element_type;
        using OIdx2 = typename decltype(idx_range_fct_parall_2)::discrete_element_type;
        OIdx1 slice_idx_1;
        OIdx2 slice_idx_2;
        /*
            If function_values is defined on the same index range than the one given
            in the collection of index range, then it means we working on the values 
            of the function.
            Otherwise,  it means we working on the first derivatives of the function. 
            The operator is then applied to compute the cross-derivatives. The index
            we return is (1) for the first derivative. 
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

        return std::make_tuple(slice_idx_1, slice_idx_2);
    }
};
