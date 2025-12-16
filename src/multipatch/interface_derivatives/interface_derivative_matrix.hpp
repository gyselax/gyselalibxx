// SPDX-License-Identifier: MIT

#pragma once

#include <typeinfo>

#include <ddc/ddc.hpp>

#include <boost/type_index.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "types.hpp"
#include "view.hpp"

template <
        class Connectivity,
        class Grid1D,
        class PatchSeq,
        ddc::BoundCond LowerBound,
        ddc::BoundCond UpperBound,
        class DerivativesCalculatorCollection>
class InterfaceExactDerivativeMatrix;

/**
  * @brief Class to compute the interface derivatives along a given direction 
  * on all the interfaces.  
  * 
  * This operator is implemented for conforming equivalent global meshes. 
  * It uses the exact formula. 
  * 
  * During the instantiation of the class, the matrix (I-M) is computed are stored in the class. 
  * When we call the operator .solve_deriv() or .solve_cross_deriv(): 
  *     - InterfaceExactDerivativeMatrix computes the vector C from the given values. 
  *     - It inverses the matrix system: (I-M)S = C. 
  *     - The interface derivatives are stored in the vector S. 
  *     - It updates the given derivatives field with the values of the vector S. 
  * 
  * @tparam Connectivity A MultipatchConnectivity class describing all the patch connections.
  * @tparam Grid1D A given direction.
  * @tparam Patches List of patches containing all the involved patches in the given direction. 
  * @tparam LowerBound Lower/left boundary condition of the first local spline along the given direction. 
  * @tparam UpperBound Upper/right boundary condition of the last local spline along the given direction. 
  * @tparam DerivativesCalculatorCollection A SingleInterfaceDerivativesCalculatorCollection
  * that stores the SingleInterfaceDerivativesCalculator needed to compute the interface derivatives. 
  */
template <
        class Connectivity,
        class Grid1D,
        class... Patches,
        ddc::BoundCond LowerBound,
        ddc::BoundCond UpperBound,
        class DerivativesCalculatorCollection>
class InterfaceExactDerivativeMatrix<
        Connectivity,
        Grid1D,
        ddc::detail::TypeSeq<Patches...>,
        LowerBound,
        UpperBound,
        DerivativesCalculatorCollection>
{
    /*
        All the interfaces given as input to the MultipatchConnectivity class.
         We expect the parameters defined on these interfaces.
    */
    using all_interface_collection = typename Connectivity::interface_collection;
    // All the sorted interfaces with the correct orientation.
    using interface_sorted_collection =
            typename Connectivity::get_all_interfaces_along_direction_t<Grid1D>;

    // All the patches of the geometry.
    using all_patches = typename Connectivity::all_patches;

    static_assert(
            is_single_derivative_calculator_collection_v<DerivativesCalculatorCollection>,
            "Please provide a SingleInterfaceDerivativesCalculatorCollection type.");

    static constexpr std::size_t number_of_interfaces
            = ddc::type_seq_size_v<interface_sorted_collection>;

    static_assert(
            ((LowerBound == ddc::BoundCond::PERIODIC) == (UpperBound == ddc::BoundCond::PERIODIC)),
            "If one boundary is periodic, the other boundary should be too.");

    static constexpr bool is_periodic = (LowerBound == ddc::BoundCond::PERIODIC);

    /*
        Remove all the interfaces with an OutsideEdge.
        Rely on get_all_interfaces_along_direction_t to order the interfaces.
    */
    using outer_interface_collection = std::conditional_t<
            is_periodic,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<
                    ddc::type_seq_element_t<0, interface_sorted_collection>,
                    ddc::type_seq_element_t<
                            number_of_interfaces - 1,
                            interface_sorted_collection>>>;

    // All the interfaces in the Grid1D direction without the interfaces with OutsideEdge.
    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_sorted_collection, outer_interface_collection>;

    // Number of inner interfaces. It fixes the number of elements in the matrix and vectors.
    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;

    using FirstInterface = ddc::type_seq_element_t<0, inner_interface_collection>;

    // First patch of the sorted interfaces.
    using FirstPatch = std::conditional_t<
            is_periodic,
            find_patch_t<Grid1D, all_patches>,
            typename FirstInterface::Edge2::associated_patch>;

    // Sequence of 1D grids in the direction of Grid1D over the patches.
    using Grid1DSeq = collect_grids_on_dim_t<
            find_patch_t<Grid1D, all_patches>,
            Grid1D,
            interface_sorted_collection>;

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

    // Tag to identify we are computing the first derivatives.
    struct eval_deriv
    {
    };

    // Tag to identify we are computing the cross-derivatives.
    struct eval_cross_deriv
    {
    };

    using Matrix = gko::matrix::Dense<double>;

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using first_patch = FirstPatch;

private:
    // Matrix (I-M) of the documentation.
    std::shared_ptr<Matrix> m_matrix;
    // Vector C of the documentation.
    std::shared_ptr<Matrix> m_vector;
    // Vector S of the documentation.
    std::shared_ptr<Matrix> m_interface_derivatives;

    // Use a conjugate gradient (CG) solver.
    using SolverCG = gko::solver::Cg<>;
    // Use a Jacobi preconditioner.
    using PreconditionerJ = gko::preconditioner::Jacobi<>;
    std::unique_ptr<SolverCG> m_solver;

    MultipatchType<IdxRangeOnPatch, Patches...> const& m_idx_ranges;

    DerivativesCalculatorCollection const& m_derivatives_calculators;

    // Vector containing the signs of the derivatives along the parallel direction of the interfaces.
    std::array<int, ddc::type_seq_size_v<Grid1DSeq>> m_vector_par_signs;

public:
    ~InterfaceExactDerivativeMatrix() {}

    /**
     * @brief Instantiate InterfaceExactDerivativeMatrix. 
     * 
     * It creates a dense Gingko matrix, a vector for the rhs and a vector for the solution. 
     * The matrix is initialised with the derivatives calculator collection given in input. 
     * 
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param derivatives_calculators SingleInterfaceDerivativesCalculatorCollection containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     * @param reduction_factor Reduction factor for a Gingko dense matrix with a conjugate gradient solver 
     *          and a Jacobi preconditioner. 
     */
    InterfaceExactDerivativeMatrix(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            DerivativesCalculatorCollection const& derivatives_calculators,
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

        set_sign_vector(
                std::make_integer_sequence<std::size_t, ddc::type_seq_size_v<Grid1DSeq> - 1> {});
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
    void solve_deriv(
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
     * @warning The global mesh is supposed to be conforming. 
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
     * @brief Solve the matrix system MS = C to determine all the interface cross-derivatives in the Grid1D
     * direction at a given index.
     * 
     * It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces. 
     *
     * @tparam IdxPar Type for the given index. It can be on the first or the second grid of the first patch.
     *
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @param[in] idx_par Index of the line in the Grid1D direction where we compute all the interface cross-derivatives.
     */
    template <class IdxPar>
    void solve_cross_deriv(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        solve<eval_cross_deriv>(functions_and_derivs, idx_par);
    }

    /**
     * @brief Solve the matrix system MS = C to determine all the corner cross-derivatives in the Grid1D
     * direction at all the last and the first indices of the first patch.
     * 
     * It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces. 
     *
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @warning The global mesh is supposed to be conforming. 
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
                "InterfaceExactDerivativeMatrix<...>::first_patch.");

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

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_same_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

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

        // Add coefficients for the periodic case.
        if constexpr (is_periodic && I == 0) {
            m_matrix->at(I, n_inner_interfaces - 1) = coefs[0];
        }
        if constexpr (is_periodic && I == n_inner_interfaces - 1) {
            m_matrix->at(I, 0) = coefs[2];
        }
    }

    /**
     * @brief Set a vector with orientation sign for the parallel grids.
     * In the computation of the cross-derivatives, we need to know if the direction 
     * of the parallel grids agree or not. It is precomputed here and stored in a vector.
     */
    template <std::size_t... I>
    void set_sign_vector(std::integer_sequence<std::size_t, I...>)
    {
        // Arbitrary default value for the first patch.
        m_vector_par_signs[0] = 1;
        (set_sign_vector_line<I>(), ...);
    }

    /// @brief Fill in the Ith line of the vector contening the orientation signs of the parallel grids.
    template <std::size_t I>
    void set_sign_vector_line()
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_same_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        constexpr Extremity extremity_1 = EquivalentInterfaceI::Edge1::extremity;
        constexpr Extremity extremity_2 = EquivalentInterfaceI::Edge2::extremity;

        using Patch_1 = typename EquivalentInterfaceI::Edge1::associated_patch;
        using Patch_2 = typename EquivalentInterfaceI::Edge2::associated_patch;

        using GridPerp1 = typename EquivalentInterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename EquivalentInterfaceI::Edge2::perpendicular_grid;

        // If the perpendicular grids are agree on one dimension (first or second).
        constexpr bool is_same_dimension
                = ((std::is_same_v<GridPerp1, typename Patch_1::Grid1>)
                   == (std::is_same_v<GridPerp2, typename Patch_2::Grid1>));

        // If the perpendicular grids agree on the direction.
        constexpr bool are_same_direction_perp = (extremity_1 != extremity_2);

        // If we need to change the sign compared to the previous interface.
        constexpr bool change_sign = (are_same_direction_perp == is_same_dimension);

        const int previous_sign = m_vector_par_signs[I];
        constexpr int changing_sign = change_sign - !change_sign;

        m_vector_par_signs[I + 1] = previous_sign * changing_sign;
    }

    /// @brief Set the vector C.
    template <typename eval_type, class GridPar1D, std::size_t... I>
    void set_vector(
            Idx<GridPar1D> const& slice_idx,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            std::integer_sequence<std::size_t, I...>)
    {
        // Check we start the loop with the first patch.
        static_assert(
                (std::is_same_v<GridPar1D, typename FirstPatch::Grid1>)
                || (std::is_same_v<GridPar1D, typename FirstPatch::Grid2>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim1>>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim2>>));

        // Locally transform the index into an integer to avoid type issues.
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
            int& slice_previous_idx_value,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_same_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        using Patch_1 = typename EquivalentInterfaceI::Edge1::associated_patch;
        using Patch_2 = typename EquivalentInterfaceI::Edge2::associated_patch;

        constexpr Extremity extremity_1 = EquivalentInterfaceI::Edge1::extremity;
        constexpr Extremity extremity_2 = EquivalentInterfaceI::Edge2::extremity;

        using GridPerp1 = typename EquivalentInterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename EquivalentInterfaceI::Edge2::perpendicular_grid;

        using DerivPerp1 = typename ddc::Deriv<typename GridPerp1::continuous_dimension_type>;
        using DerivPerp2 = typename ddc::Deriv<typename GridPerp2::continuous_dimension_type>;

        // Type for the slice_previous_idx_value index.
        using PatchSliceIdx = std::conditional_t<is_same_orientation, Patch_1, Patch_2>;
        auto [idx_slice_1, idx_slice_2]
                = get_slice_indexes<EquivalentInterfaceI, PatchSliceIdx>(slice_previous_idx_value);

        DerivFieldOnPatch_host<Patch_1> function_and_derivs_1
                = functions_and_derivs.template get<Patch_1>();
        DerivFieldOnPatch_host<Patch_2> function_and_derivs_2
                = functions_and_derivs.template get<Patch_2>();

        // Use the function values to compute the first derivatives.
        DField<typename Patch_1::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_1
                = function_and_derivs_1.get_values_field();
        DField<typename Patch_2::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_2
                = function_and_derivs_2.get_values_field();

        // Compute the coefficient c_I for the interface I.
        double const lin_comb_funct
                = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                          .get_function_coefficients(
                                  get_const_field(function_and_derivs_1[idx_slice_1]),
                                  get_const_field(function_and_derivs_2[idx_slice_2]));

        /*
            If the orientation of the interface is not correct, we change the sign of the sum.
            Because c_{I 1|2} = - c_{I 2|1}.
        */
        constexpr int sign = is_same_orientation - !is_same_orientation;

        m_vector->get_values()[I] = sign * lin_comb_funct;

        // Add the boundary derivatives for global Hermite boundary conditions. ------------------
        constexpr bool is_lower_bound_deriv_dependent
                = (LowerBound == ddc::BoundCond::HERMITE && I == 0);
        constexpr bool is_upper_bound_deriv_dependent
                = (UpperBound == ddc::BoundCond::HERMITE && I == n_inner_interfaces - 1);

        if constexpr (is_lower_bound_deriv_dependent || is_upper_bound_deriv_dependent) {
            IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch_1>());
            IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch_2>());

            IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1
                    = function_and_derivs_1.template idx_range_for_deriv<GridPerp1>();
            IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2
                    = function_and_derivs_2.template idx_range_for_deriv<GridPerp2>();

            Idx<GridPerp1> idx_deriv_1
                    = get_idx_other_interface(idx_range_slice_dperp_1, extremity_1);
            Idx<GridPerp2> idx_deriv_2
                    = get_idx_other_interface(idx_range_slice_dperp_2, extremity_2);

            Idx<DerivPerp1, GridPerp1> idx_slice_deriv_1(Idx<DerivPerp1>(1), idx_deriv_1);
            Idx<DerivPerp2, GridPerp2> idx_slice_deriv_2(Idx<DerivPerp2>(1), idx_deriv_2);

            // If the directions of the two perpendicular grids desagree, we will change the sign.
            constexpr bool are_same_direction_perp = (extremity_1 != extremity_2);
            constexpr int sign_deriv_perp_2 = are_same_direction_perp - !are_same_direction_perp;

            // If Patch1 follows the orientation of the sorted 1D grid sequence.
            constexpr bool is_per_1_well_oriented
                    = (is_same_orientation && (extremity_1 == Extremity::BACK))
                      || (!is_same_orientation && (extremity_1 == Extremity::FRONT));
            constexpr int sign_deriv_perp_1 = is_per_1_well_oriented - !is_per_1_well_oriented;

            // Change the sign of the derivatives if Patch1 is ill-oriented.
            const double deriv_1
                    = function_and_derivs_1[idx_slice_deriv_1](idx_slice_1) * sign_deriv_perp_1;
            // Change in addition the sign of the derivative on Patch2, if the grid directions desagree.
            const double deriv_2 = function_and_derivs_2[idx_slice_deriv_2](idx_slice_2)
                                   * sign_deriv_perp_1 * sign_deriv_perp_2;

            const double coeff_deriv_1
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_1();
            const double coeff_deriv_2
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_2();

            if constexpr (is_lower_bound_deriv_dependent) {
                m_vector->get_values()[I]
                        += is_same_orientation ? coeff_deriv_1 * deriv_1 : coeff_deriv_2 * deriv_2;
            } else {
                m_vector->get_values()[I]
                        += is_same_orientation ? coeff_deriv_2 * deriv_2 : coeff_deriv_1 * deriv_1;
            }
        }
    }

    /// @brief Fill in the Ith line of the vector C corresponding.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_cross_deriv>, bool> = true>
    void set_line_vector(
            int& slice_previous_idx_value,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_same_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        constexpr Extremity extremity_1 = EquivalentInterfaceI::Edge1::extremity;
        constexpr Extremity extremity_2 = EquivalentInterfaceI::Edge2::extremity;

        using Patch_1 = typename EquivalentInterfaceI::Edge1::associated_patch;
        using Patch_2 = typename EquivalentInterfaceI::Edge2::associated_patch;

        // Type for the slice_previous_idx_value index.
        using PatchSliceIdx = std::conditional_t<is_same_orientation, Patch_1, Patch_2>;
        auto [idx_slice_1, idx_slice_2]
                = get_slice_indexes<EquivalentInterfaceI, PatchSliceIdx>(slice_previous_idx_value);

        DerivFieldOnPatch_host<Patch_1> function_and_derivs_1
                = functions_and_derivs.template get<Patch_1>();
        DerivFieldOnPatch_host<Patch_2> function_and_derivs_2
                = functions_and_derivs.template get<Patch_2>();

        // Use the first derivatives to compute the cross-derivatives.
        using GridPerp1 = typename EquivalentInterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename EquivalentInterfaceI::Edge2::perpendicular_grid;

        using GridPar1 = typename EquivalentInterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename EquivalentInterfaceI::Edge2::parallel_grid;

        using DerivPar1 = typename ddc::Deriv<typename GridPar1::continuous_dimension_type>;
        using DerivPar2 = typename ddc::Deriv<typename GridPar2::continuous_dimension_type>;

        IdxRange<GridPar1> idx_range_par_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<GridPar2> idx_range_par_2(m_idx_ranges.template get<Patch_2>());

        // The slice indices has to be a point at a corner, i.e. index range boundaries.
        assert((idx_slice_1 == idx_range_par_1.front()) || (idx_slice_1 == idx_range_par_1.back()));
        assert((idx_slice_2 == idx_range_par_2.front()) || (idx_slice_2 == idx_range_par_2.back()));

        IdxRangeSlice<GridPar1> idx_range_slice_dpar_1
                = function_and_derivs_1.template idx_range_for_deriv<GridPar1>();
        IdxRangeSlice<GridPar2> idx_range_slice_dpar_2
                = function_and_derivs_2.template idx_range_for_deriv<GridPar2>();

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

        // Compute the coefficient c_I for the interface I.
        double const lin_comb_funct = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                                              .get_derivatives_coefficients(
                                                      get_const_field(derivs_1),
                                                      get_const_field(derivs_2));

        // If the orientations are not the same, we change the sign of the sum.
        const int sign = is_same_orientation - !is_same_orientation;
        const int sign_deriv_par_1 = m_vector_par_signs[ddc::type_seq_rank_v<GridPerp1, Grid1DSeq>];

        /*
            We need to change the sign if the orientation of the interface desagrees with 
            the sorted 1D grid sequence. Because,
                    c_{I 1|2} = - c_{I 2|1}.
            In get_derivatives_coefficients(), the sign of the derivatives on Patch2 is changed
            if the edge orientation desagree. We need to correct the sign, if Patch1 was 
            ill-oriented. Because, 
                    d(f(-x)) = - df(-x).
        */
        m_vector->get_values()[I] = sign * sign_deriv_par_1 * lin_comb_funct;

        // Add the boundary derivatives for global Hermite boundary conditions. ------------------
        constexpr bool is_lower_bound_deriv_dependent
                = (LowerBound == ddc::BoundCond::HERMITE && I == 0);
        constexpr bool is_upper_bound_deriv_dependent
                = (UpperBound == ddc::BoundCond::HERMITE && I == n_inner_interfaces - 1);

        if constexpr (is_lower_bound_deriv_dependent || is_upper_bound_deriv_dependent) {
            using Grid1_1 = typename Patch_1::Grid1;
            using Grid2_1 = typename Patch_1::Grid2;
            using Grid1_2 = typename Patch_2::Grid1;
            using Grid2_2 = typename Patch_2::Grid2;

            using Deriv1_1 = ddc::Deriv<typename Patch_1::Dim1>;
            using Deriv2_1 = ddc::Deriv<typename Patch_1::Dim2>;
            using Deriv1_2 = ddc::Deriv<typename Patch_2::Dim1>;
            using Deriv2_2 = ddc::Deriv<typename Patch_2::Dim2>;

            IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch_1>());
            IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch_2>());

            IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1
                    = function_and_derivs_1.template idx_range_for_deriv<GridPerp1>();
            IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2
                    = function_and_derivs_2.template idx_range_for_deriv<GridPerp2>();

            Idx<GridPerp1> idx_deriv_perp_1
                    = get_idx_other_interface(idx_range_slice_dperp_1, extremity_1);
            Idx<GridPerp2> idx_deriv_perp_2
                    = get_idx_other_interface(idx_range_slice_dperp_2, extremity_2);

            Idx<Grid1_1, Grid2_1> idx_d1d2_1(idx_deriv_par_1, idx_deriv_perp_1);
            Idx<Grid1_2, Grid2_2> idx_d1d2_2(idx_deriv_par_2, idx_deriv_perp_2);

            Idx<Grid1_1> idx_d1_1(idx_d1d2_1);
            Idx<Grid2_1> idx_d2_1(idx_d1d2_1);
            Idx<Grid1_2> idx_d1_2(idx_d1d2_2);
            Idx<Grid2_2> idx_d2_2(idx_d1d2_2);

            Idx<Deriv1_1, Grid1_1, Deriv2_1, Grid2_1>
                    idx_cross_deriv1(Idx<Deriv1_1>(1), idx_d1_1, Idx<Deriv2_1>(1), idx_d2_1);
            Idx<Deriv1_2, Grid1_2, Deriv2_2, Grid2_2>
                    idx_cross_deriv2(Idx<Deriv1_2>(1), idx_d1_2, Idx<Deriv2_2>(1), idx_d2_2);

            // If Patch1 follows the orientation of the sorted 1D grid sequence.
            constexpr bool is_per_1_well_oriented
                    = (is_same_orientation && (extremity_1 == Extremity::BACK))
                      || (!is_same_orientation && (extremity_1 == Extremity::FRONT));
            constexpr int sign_deriv_perp_1 = is_per_1_well_oriented - !is_per_1_well_oriented;

            // Sign for the parallel axis.
            const int sign_deriv_par_1
                    = m_vector_par_signs[ddc::type_seq_rank_v<GridPerp1, Grid1DSeq>];
            // If the perpendicular grids of the two patches agree.
            constexpr bool are_same_direction_perp = (extremity_1 != extremity_2);
            // If the parallel grids of the two patches agree.
            constexpr bool are_same_direction_par = EquivalentInterfaceI::orientations_agree;

            constexpr int sign_deriv_perp_2 = are_same_direction_perp - !are_same_direction_perp;
            constexpr int sign_deriv_par_2 = are_same_direction_par - !are_same_direction_par;

            /*
                Change the sign of the cross-derivative if Patch1 has its perpendicular grid ill-oriented.
                Change the sign of the cross-derivative if Patch1 has its parallel grid ill-oriented. 
            */
            const double cross_deriv_1 = function_and_derivs_1(idx_cross_deriv1) * sign_deriv_perp_1
                                         * sign_deriv_par_1;
            /*
                Change the sign of the cross-derivative if Patch1 is ill-oriented.
                Correct the sign of the cross-derivative from the parallel axis of Patch1.
                Change the sign of the cross-derivative if the orientation of the perpendicular grid 
                of Patch2 desagrees with Patch1.
                Change the sign of the cross-derivative if the orientation of the parallel grid 
                of Patch2 desagrees with Patch1.
            */
            const double cross_deriv_2 = function_and_derivs_2(idx_cross_deriv2) * sign_deriv_perp_1
                                         * sign_deriv_par_1 * sign_deriv_perp_2 * sign_deriv_par_2;

            const double coeff_deriv_1
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_1();
            const double coeff_deriv_2
                    = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                              .get_coeff_deriv_patch_2();

            if constexpr (is_lower_bound_deriv_dependent) {
                m_vector->get_values()[I] += is_same_orientation ? coeff_deriv_1 * cross_deriv_1
                                                                 : coeff_deriv_2 * cross_deriv_2;
            } else {
                m_vector->get_values()[I] += is_same_orientation ? coeff_deriv_2 * cross_deriv_2
                                                                 : coeff_deriv_1 * cross_deriv_1;
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

        // Locally transform the index into an integer to avoid type issues.
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
            int& slice_previous_idx_value)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_same_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        using Patch_1 = typename EquivalentInterfaceI::Edge1::associated_patch;
        using Patch_2 = typename EquivalentInterfaceI::Edge2::associated_patch;

        constexpr Extremity extremity_1 = EquivalentInterfaceI::Edge1::extremity;
        constexpr Extremity extremity_2 = EquivalentInterfaceI::Edge2::extremity;

        using GridPerp1 = typename EquivalentInterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename EquivalentInterfaceI::Edge2::perpendicular_grid;

        using GridPar1 = typename EquivalentInterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename EquivalentInterfaceI::Edge2::parallel_grid;

        using DerivPerp1 = typename ddc::Deriv<typename GridPerp1::continuous_dimension_type>;
        using DerivPerp2 = typename ddc::Deriv<typename GridPerp2::continuous_dimension_type>;

        // Get the fields of the left and right patch of the interface.
        DerivFieldOnPatch_host<Patch_1> function_and_derivs_1
                = functions_and_derivs.template get<Patch_1>();
        DerivFieldOnPatch_host<Patch_2> function_and_derivs_2
                = functions_and_derivs.template get<Patch_2>();

        // Get the correct index ranges and indices for the slices.
        IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch_2>());

        IdxRangeSlice<GridPerp1> idx_range_slice_dperp_1
                = function_and_derivs_1.template idx_range_for_deriv<GridPerp1>();
        IdxRangeSlice<GridPerp2> idx_range_slice_dperp_2
                = function_and_derivs_2.template idx_range_for_deriv<GridPerp2>();

        Idx<GridPerp1> idx_deriv_1 = get_idx_interface(idx_range_slice_dperp_1, extremity_1);
        Idx<GridPerp2> idx_deriv_2 = get_idx_interface(idx_range_slice_dperp_2, extremity_2);

        Idx<DerivPerp1, GridPerp1> idx_slice_deriv_1(Idx<DerivPerp1>(1), idx_deriv_1);
        Idx<DerivPerp2, GridPerp2> idx_slice_deriv_2(Idx<DerivPerp2>(1), idx_deriv_2);

        // Type for the slice_previous_idx_value index.
        using PatchSliceIdx = std::conditional_t<is_same_orientation, Patch_1, Patch_2>;
        auto [idx_slice_1, idx_slice_2]
                = get_slice_indexes<EquivalentInterfaceI, PatchSliceIdx>(slice_previous_idx_value);

        // Select the correct data to update.
        DField<IdxRange<GridPar1>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_1
                = function_and_derivs_1[idx_slice_deriv_1];
        DField<IdxRange<GridPar2>, Kokkos::HostSpace, Kokkos::layout_stride> deriv_2
                = function_and_derivs_2[idx_slice_deriv_2];

        // The orientation of Patch1 follows the global orientation on the sorted 1D grid sequence.
        constexpr bool is_per_1_well_oriented
                = (is_same_orientation && (extremity_1 == Extremity::BACK))
                  || (!is_same_orientation && (extremity_1 == Extremity::FRONT));
        constexpr int sign_deriv_perp_1 = is_per_1_well_oriented - !is_per_1_well_oriented;

        // The direction of the two perpendicular grids of the patches agree.
        constexpr bool are_patches_same_direction = (extremity_1 != extremity_2);
        constexpr int sign_change_direction
                = are_patches_same_direction - !are_patches_same_direction;

        // Change again the sign of the derivatives if Patch1 was ill-oriented.
        deriv_1(idx_slice_1) = m_interface_derivatives->get_values()[I] * sign_deriv_perp_1;
        // Change again the sign of the derivative on Patch2, if the grid directions desagreed.
        deriv_2(idx_slice_2) = m_interface_derivatives->get_values()[I] * sign_deriv_perp_1
                               * sign_change_direction;
    }

    /// @brief Associate the Ith derivative values to the correct cross-derivative.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_cross_deriv>, bool> = true>
    void update_derivatives_at_interface(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& slice_previous_idx_value)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_same_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        using Patch_1 = typename EquivalentInterfaceI::Edge1::associated_patch;
        using Patch_2 = typename EquivalentInterfaceI::Edge2::associated_patch;

        constexpr Extremity extremity_1 = EquivalentInterfaceI::Edge1::extremity;
        constexpr Extremity extremity_2 = EquivalentInterfaceI::Edge2::extremity;

        using GridPerp1 = typename EquivalentInterfaceI::Edge1::perpendicular_grid;

        using GridPar1 = typename EquivalentInterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename EquivalentInterfaceI::Edge2::parallel_grid;

        using Grid1_1 = typename Patch_1::Grid1;
        using Grid2_1 = typename Patch_1::Grid2;
        using Grid1_2 = typename Patch_2::Grid1;
        using Grid2_2 = typename Patch_2::Grid2;

        using Deriv1_1 = ddc::Deriv<typename Patch_1::Dim1>;
        using Deriv2_1 = ddc::Deriv<typename Patch_1::Dim2>;
        using Deriv1_2 = ddc::Deriv<typename Patch_2::Dim1>;
        using Deriv2_2 = ddc::Deriv<typename Patch_2::Dim2>;

        constexpr bool is_grid_par_1_on_dim1 = std::is_same_v<GridPar1, Grid1_1>;
        constexpr bool is_grid_par_2_on_dim1 = std::is_same_v<GridPar2, Grid1_2>;

        // Get the fields of the left and right patch of the interface.
        DerivFieldOnPatch_host<Patch_1> function_and_derivs_1
                = functions_and_derivs.template get<Patch_1>();
        DerivFieldOnPatch_host<Patch_2> function_and_derivs_2
                = functions_and_derivs.template get<Patch_2>();

        // Get the correct index ranges and indices for the slices.
        IdxRange<Grid1_1> idx_range_grid1_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<Grid2_1> idx_range_grid2_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<Grid1_2> idx_range_grid1_2(m_idx_ranges.template get<Patch_2>());
        IdxRange<Grid2_2> idx_range_grid2_2(m_idx_ranges.template get<Patch_2>());

        IdxRangeSlice<Grid1_1> idx_range_slice_d1_1
                = function_and_derivs_1.template idx_range_for_deriv<Grid1_1>();
        IdxRangeSlice<Grid2_1> idx_range_slice_d2_1
                = function_and_derivs_1.template idx_range_for_deriv<Grid2_1>();
        IdxRangeSlice<Grid1_2> idx_range_slice_d1_2
                = function_and_derivs_2.template idx_range_for_deriv<Grid1_2>();
        IdxRangeSlice<Grid2_2> idx_range_slice_d2_2
                = function_and_derivs_2.template idx_range_for_deriv<Grid2_2>();

        IdxRange<GridPar1> idx_range_par_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<GridPar2> idx_range_par_2(m_idx_ranges.template get<Patch_2>());

        // Type for the slice_previous_idx_value index.
        using PatchSliceIdx = std::conditional_t<is_same_orientation, Patch_1, Patch_2>;
        auto [idx_slice_1, idx_slice_2]
                = get_slice_indexes<EquivalentInterfaceI, PatchSliceIdx>(slice_previous_idx_value);

        // The slice indices has to be a point at a corner, i.e. index range boundaries.
        assert((idx_slice_1 == idx_range_par_1.front()) || (idx_slice_1 == idx_range_par_1.back()));
        assert((idx_slice_2 == idx_range_par_2.front()) || (idx_slice_2 == idx_range_par_2.back()));

        const bool is_idx_par_min_1 = (idx_slice_1 == idx_range_par_1.front());
        const bool is_idx_par_min_2 = (idx_slice_2 == idx_range_par_2.front());

        // If Patch1 follows the orientation of the sorted 1D grid sequence.
        constexpr bool is_per_1_well_oriented
                = (is_same_orientation && (extremity_1 == Extremity::BACK))
                  || (!is_same_orientation && (extremity_1 == Extremity::FRONT));
        constexpr int sign_deriv_perp_1 = is_per_1_well_oriented - !is_per_1_well_oriented;

        // Sign for the parallel axis.
        const int sign_deriv_par_1 = m_vector_par_signs[ddc::type_seq_rank_v<GridPerp1, Grid1DSeq>];

        // If the perpendicular grids of the two patches agree.
        constexpr bool are_same_direction_perp = (extremity_1 != extremity_2);
        // If the parallel grids of the two patches agree.
        constexpr bool are_same_direction_par = EquivalentInterfaceI::orientations_agree;

        constexpr int sign_deriv_perp_2 = are_same_direction_perp - !are_same_direction_perp;
        constexpr int sign_deriv_par_2 = are_same_direction_par - !are_same_direction_par;

        // --- Update the cross-derivative on Patch_1
        /*
            Change again the sign of the cross-derivative if Patch1 had its perpendicular grid ill-oriented.
            Change again the sign of the cross-derivative if Patch1 had its parallel grid ill-oriented. 
        */
        Idx<Deriv1_1, Grid1_1, Deriv2_1, Grid2_1> idx_cross_deriv1
                = get_corner_idx<Deriv1_1, Grid1_1, Deriv2_1, Grid2_1>(
                        extremity_1,
                        is_grid_par_1_on_dim1,
                        is_idx_par_min_1,
                        idx_range_slice_d1_1,
                        idx_range_slice_d2_1);
        function_and_derivs_1(idx_cross_deriv1)
                = m_interface_derivatives->get_values()[I] * sign_deriv_perp_1 * sign_deriv_par_1;

        // --- Update the cross-derivative on Patch_2
        /*
            Change again the sign of the cross-derivative if Patch1 was ill-oriented.
            Correct again the sign of the cross-derivative from the parallel axis of Patch1.
            Change again the sign of the cross-derivative if the orientation of the perpendicular grid 
            of Patch2 desagreed with Patch1.
            Change again the sign of the cross-derivative if the orientation of the parallel grid 
            of Patch2 desagreed with Patch1.
        */
        Idx<Deriv1_2, Grid1_2, Deriv2_2, Grid2_2> idx_cross_deriv2
                = get_corner_idx<Deriv1_2, Grid1_2, Deriv2_2, Grid2_2>(
                        extremity_2,
                        is_grid_par_2_on_dim1,
                        is_idx_par_min_2,
                        idx_range_slice_d1_2,
                        idx_range_slice_d2_2);
        function_and_derivs_2(idx_cross_deriv2) = m_interface_derivatives->get_values()[I]
                                                  * sign_deriv_perp_1 * sign_deriv_par_1
                                                  * sign_deriv_perp_2 * sign_deriv_par_2;
    }


    /**
     * @brief Get the indices on Patch1 and Patch2 of the given interface.
     * The index on Patch1 corresponds to the value given in input
     * (if the interface is well oriented, otherwise it is Patch2).
     * The index on Patch2 is the equivalent index of the one on Patch1.
     * It also updates the value given in input with the value of the index
     * on Patch2 for the next interface 
     * (if the next interface is well oriented, otherwise it is Patch1).
     * 
     * Example for well oriented interfaces:
     * 
     *  |           Patch1 | Patch2      Patch1 | Patch2        |  
     *  |            idx_1 > idx_2        idx_1 > idx_2         |
     *  |              ^   |   v            ^   |   v           |
     *  |idx_val   idx_val | idx_val    idx_val | idx_val       |
     */
    template <typename InterfaceI, class PatchSliceIdx>
    auto get_slice_indexes(int& slice_idx_value)
    {
        using Patch_1 = typename InterfaceI::Edge1::associated_patch;
        using Patch_2 = typename InterfaceI::Edge2::associated_patch;

        using ParallGrid1 = typename InterfaceI::Edge1::parallel_grid;
        using ParallGrid2 = typename InterfaceI::Edge2::parallel_grid;

        IdxRange<ParallGrid1> const idx_range_parall_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<ParallGrid2> const idx_range_parall_2(m_idx_ranges.template get<Patch_2>());

        EdgeTransformation<InterfaceI> index_converter(idx_range_parall_1, idx_range_parall_2);

        if constexpr (std::is_same_v<Patch_1, PatchSliceIdx>) {
            Idx<ParallGrid1> slice_idx_1(slice_idx_value);
            Idx<ParallGrid2> slice_idx_2(index_converter(slice_idx_1));
            slice_idx_value = (slice_idx_2 - idx_range_parall_2.front()).value();
            return std::make_tuple(slice_idx_1, slice_idx_2);
        } else {
            Idx<ParallGrid2> slice_idx_2(slice_idx_value);
            Idx<ParallGrid1> slice_idx_1(index_converter(slice_idx_2));
            slice_idx_value = (slice_idx_1 - idx_range_parall_1.front()).value();
            return std::make_tuple(slice_idx_1, slice_idx_2);
        }
    }

    /**
     * @brief Get the index of the given index range slice corresponding to the 
     * FRONT extremity of the 1D grid. 
     */
    template <class Grid1Or2>
    Idx<Grid1Or2> get_idx_interface(
            IdxRangeSlice<Grid1Or2> const idx_range_slice,
            Extremity extremity) const
    {
        return (extremity == Extremity::FRONT) ? idx_range_slice.front() : idx_range_slice.back();
    }

    /**
     * @brief Get the index of the given index range slice corresponding to the 
     * BACK extremity of the 1D grid. 
     */
    template <class Grid1Or2>
    Idx<Grid1Or2> get_idx_other_interface(
            IdxRangeSlice<Grid1Or2> const idx_range_slice,
            Extremity extremity) const
    {
        return (extremity == Extremity::BACK) ? idx_range_slice.front() : idx_range_slice.back();
    }

    /**
     * @brief Get the index to extract the cross-derivative corresponding to the 
     * correct corner of the patch. 
     */
    template <class Deriv1, class Grid1, class Deriv2, class Grid2>
    Idx<Deriv1, Grid1, Deriv2, Grid2> get_corner_idx(
            Extremity const extremity,
            bool const is_grid_par_on_dim1,
            bool const is_idx_par_min,
            IdxRangeSlice<Grid1> const idx_range_slice_d1,
            IdxRangeSlice<Grid2> const idx_range_slice_d2) const
    {
        Idx<Grid1> idx_d1;
        Idx<Grid2> idx_d2;
        if (is_idx_par_min && (extremity == Extremity::FRONT)) {
            // corner min min
            idx_d1 = idx_range_slice_d1.front();
            idx_d2 = idx_range_slice_d2.front();
        } else if (!is_idx_par_min && (extremity == Extremity::BACK)) {
            // corner max max
            idx_d1 = idx_range_slice_d1.back();
            idx_d2 = idx_range_slice_d2.back();
        } else if (
                (is_grid_par_on_dim1 && (extremity == Extremity::BACK))
                || (!is_grid_par_on_dim1 && (extremity == Extremity::FRONT))) {
            // corner min max
            idx_d1 = idx_range_slice_d1.front();
            idx_d2 = idx_range_slice_d2.back();
        } else {
            // corner max min
            idx_d1 = idx_range_slice_d1.back();
            idx_d2 = idx_range_slice_d2.front();
        }
        return Idx<Deriv1, Grid1, Deriv2, Grid2>(Idx<Deriv1>(1), idx_d1, Idx<Deriv2>(1), idx_d2);
    }
};
