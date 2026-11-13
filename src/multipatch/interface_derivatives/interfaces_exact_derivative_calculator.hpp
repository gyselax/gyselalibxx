// SPDX-License-Identifier: MIT

#pragma once

#include <typeinfo>
#include <utility>
#include <vector>

#include <ddc/ddc.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "derivative_calculator_common.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "interface_derivative_coefficients_collection.hpp"
#include "matrix_dense.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"
#include "view.hpp"

template <
        class Connectivity,
        class Grid1D,
        class PatchSeq,
        ddc::SplineBuilderClosure LowerBound,
        ddc::SplineBuilderClosure UpperBound,
        class DerivCoeffCollection>
class InterfacesExactDerivativeCalculator;

/**
  * @brief Class to compute the interface derivatives along a given direction 
  * on all the interfaces.  
  * 
  * This operator is implemented for conforming equivalent global meshes. 
  * It uses the exact formula. 
  * 
  * During the instantiation of the class, the matrix (I-M) is computed are stored in the class. 
  * When we call operator(): 
  *     - InterfacesExactDerivativeCalculator computes the vector C from the given values. 
  *     - It inverses the matrix system: (I-M)S = C. 
  *     - The interface derivatives are stored in the vector S. 
  *     - It updates the given derivatives field with the values of the vector S. 
  * 
  * @tparam Connectivity A MultipatchConnectivity class describing all the patch connections.
  * @tparam Grid1D A given direction.
  * @tparam Patches List of patches containing all the involved patches in the given direction. 
  * @tparam LowerBound Lower/left boundary condition of the first local spline along the given direction. 
  * @tparam UpperBound Upper/right boundary condition of the last local spline along the given direction. 
  * @tparam DerivCoeffCollection A InterfaceDerivCoeffsCollection
  * that stores the SingleInterfaceDerivativesCalculator needed to compute the interface derivatives. 
  */
template <
        class Connectivity,
        class Grid1D,
        class... Patches,
        ddc::SplineBuilderClosure LowerBound,
        ddc::SplineBuilderClosure UpperBound,
        class DerivCoeffCollection>
class InterfacesExactDerivativeCalculator<
        Connectivity,
        Grid1D,
        ddc::TypeSeq<Patches...>,
        LowerBound,
        UpperBound,
        DerivCoeffCollection>
    : public DerivativeCalculatorCommon<
              Connectivity,
              Grid1D,
              ddc::TypeSeq<Patches...>,
              DerivCoeffCollection>
{
    using base_type = DerivativeCalculatorCommon<
            Connectivity,
            Grid1D,
            ddc::TypeSeq<Patches...>,
            DerivCoeffCollection>;

private:
    template <std::size_t I>
    using stored_interface = typename base_type::template stored_interface<I>;

    // The order of Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
    template <std::size_t I>
    static constexpr bool patch_order_matches_interfaces
            = base_type::template patch_order_matches_interfaces<I>;

    template <std::size_t I>
    using patch1 = typename base_type::template patch1<I>;
    template <std::size_t I>
    using patch2 = typename base_type::template patch2<I>;

    template <std::size_t I>
    using grid_perp1 = typename base_type::template grid_perp1<I>;
    template <std::size_t I>
    using grid_perp2 = typename base_type::template grid_perp2<I>;

    template <std::size_t I>
    using grid_par1 = typename base_type::template grid_par1<I>;
    template <std::size_t I>
    using grid_par2 = typename base_type::template grid_par2<I>;

    template <std::size_t I>
    using deriv_perp1 = typename base_type::template deriv_perp1<I>;
    template <std::size_t I>
    using deriv_perp2 = typename base_type::template deriv_perp2<I>;

    template <std::size_t I>
    using deriv_par1 = typename base_type::template deriv_par1<I>;
    template <std::size_t I>
    using deriv_par2 = typename base_type::template deriv_par2<I>;

    template <std::size_t I>
    static constexpr Extremity extremity_1 = base_type::template extremity_1<I>;
    template <std::size_t I>
    static constexpr Extremity extremity_2 = base_type::template extremity_2<I>;

    using base_type::n_inner_interfaces;
    using base_type::to_sign;

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using typename base_type::patch_first;

    /// @brief Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed.
    using grid_direction = Grid1D;

    /// @brief Patch of the given direction grid.
    using typename base_type::patch_direction;

    /// @brief Derivative dimension of the direction patch perpendicular to the mesh line.
    using typename base_type::deriv_perp_dim;
    /// @brief Derivative dimension of the direction patch parallel to the mesh line.
    using typename base_type::deriv_par_dim;

    /// @brief Index on the derivative dimension perpendicular to the mesh line.
    using typename base_type::idx_deriv_perp;
    /// @brief Index on the two derivative dimensions of patch_direction.
    using typename base_type::idx_cross_deriv;

private:
    using base_type::m_deriv_coeffs_collect;
    using base_type::m_idx_ranges;

    // Matrix (I-M) of the documentation.
    Matrix_Dense m_matrix;
    // Vector C of the documentation.
    std::vector<double> m_vector;
    // Vector S of the documentation.
    std::vector<double> m_interface_derivatives;

    // Vector containing the signs of the derivatives along the parallel direction of the interfaces.
    std::array<int, ddc::type_seq_size_v<typename base_type::Grid1DSeq>> m_vector_par_signs;

public:
    /**
     * @brief Instantiate InterfacesExactDerivativeCalculator. 
     * 
     * It creates a dense matrix, a vector for the rhs and a vector for the solution.
     * The matrix is initialised with the derivatives calculator collection given in input
     * and factorised with a LU decomposition.
     * 
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param deriv_coeffs_collect SingleInterfaceDerivCoeffCollection containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     */
    InterfacesExactDerivativeCalculator(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            DerivCoeffCollection const& deriv_coeffs_collect)
        : base_type(idx_ranges, deriv_coeffs_collect)
        , m_matrix(n_inner_interfaces)
        , m_vector(n_inner_interfaces)
        , m_interface_derivatives(n_inner_interfaces)
    {
        set_matrix(std::make_index_sequence<n_inner_interfaces> {});

        // Factorise the matrix to prepare the resolutions of the systems (I-M)S = C.
        m_matrix.factorise();

        set_sign_vector(std::make_index_sequence<
                        ddc::type_seq_size_v<typename base_type::Grid1DSeq> - 1> {});
    }

    /**
     * @brief Solve the matrix system (I-M)S = C to determine all the interface derivatives or cross-derivatives
     * in the Grid1D direction at a given index.
     * 
     * For the derivatives, it uses the function values to compute the first-derivatives perpendicular to the interfaces. 
     * For the cross-derivatives, it uses the first-derivatives to compute the cross-derivatives 
     * in a perpendicular direction to the interfaces. 
     *
     * @tparam IdxDeriv Type for the given derivative order, either on the perpendicular direction for the 
     *      first derivatives or both for the cross-derivatives. 
     * @tparam IdxPar Type for the given index. It can be on the first or the second grid of the first patch.
     * @param[in] deriv_order Indicate if we compute the first derivatives or the cross-derivatives. 
     *      Its value is not taken into account. 
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @param[in] idx_par Index of the line in the Grid1D direction where we compute all the interface derivatives.
     */
    template <typename IdxDeriv, class IdxPar>
    void operator()(
            IdxDeriv const& deriv_order,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        // Check the input types.
        static_assert(
                (std::is_same_v<IdxPar, typename patch_first::Idx1>)
                        || (std::is_same_v<IdxPar, typename patch_first::Idx2>),
                "The given index has to be a 1D index defined on the first patch. See "
                "InterfacesExactDerivativeCalculator<...>::patch_first.");

        // Update m_vector.
        set_vector(
                deriv_order,
                idx_par,
                functions_and_derivs,
                std::make_index_sequence<n_inner_interfaces> {});

        // Update m_interface_derivatives by solving the matrix system (I-M)S = C.
        m_interface_derivatives = m_vector;
        m_matrix.solve_inplace(
                DSpan1D(m_interface_derivatives.data(), m_interface_derivatives.size()));

        // Update derivatives.
        update_derivatives(
                deriv_order,
                functions_and_derivs,
                idx_par,
                std::make_index_sequence<n_inner_interfaces> {});
    }

    /**
     * @brief Solve the matrix system MS = C to determine all the interface derivatives in the Grid1D
     * direction at all the indices of the first patch, or all the corner cross-derivatives in the Grid1D
     * direction at all the last and the first indices of the first patch.
     * 
     * For the derivatives, it uses the function values to compute the first-derivatives perpendicular to the interfaces. 
     * For the cross-derivatives, it uses the first-derivatives to compute the cross-derivatives 
     * in a perpendicular direction to the interfaces. 
     * 
     * @tparam IdxDeriv Type for the given derivative order, either on the perpendicular direction for the 
     *      first derivatives or both for the cross-derivatives. 
     * @param[in] deriv_order Indicate if we compute the first derivatives or the cross-derivatives. 
     *      Its value is not taken into account. 
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives 
     *      we want to compute.
     * @warning The global mesh is supposed to be conforming. 
     */
    template <class IdxDeriv>
    void operator()(
            IdxDeriv const& deriv_order,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        static_assert(
                (std::is_same_v<IdxDeriv, idx_deriv_perp>)
                        || (std::is_same_v<IdxDeriv, Idx<deriv_perp_dim, deriv_par_dim>>)
                        || (std::is_same_v<IdxDeriv, Idx<deriv_par_dim, deriv_perp_dim>>),
                "The derivative order has to be on "
                "InterfacesExactDerivativeCalculator<...>::deriv_perp_dim "
                "or on InterfacesExactDerivativeCalculator<...>::deriv_perp_dim "
                "and InterfacesExactDerivativeCalculator<...>::deriv_par_dim.");

        constexpr bool is_first_idx_range_perp_on_dim1
                = ddc::in_tags_v<typename patch_first::Grid1, typename base_type::Grid1DSeq>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename patch_first::IdxRange2,
                typename patch_first::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first(m_idx_ranges.template get<patch_first>());

        if constexpr (std::is_same_v<IdxDeriv, idx_deriv_perp>) {
            ddc::host_for_each(
                    idx_range_par_first,
                    [&](typename IdxRangeParFirstType::discrete_element_type const& idx_par) {
                        (*this)(deriv_order, functions_and_derivs, idx_par);
                    });
        } else {
            (*this)(deriv_order, functions_and_derivs, idx_range_par_first.front());
            (*this)(deriv_order, functions_and_derivs, idx_range_par_first.back());
        }
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
        double const coeff_deriv_patch1 = m_deriv_coeffs_collect.template get<stored_interface<I>>()
                                                  .get_coeff_deriv_patch_1();
        double const coeff_deriv_patch2 = m_deriv_coeffs_collect.template get<stored_interface<I>>()
                                                  .get_coeff_deriv_patch_2();
        const double coeff_left
                = patch_order_matches_interfaces<I> ? coeff_deriv_patch1 : coeff_deriv_patch2;
        const double coeff_right
                = patch_order_matches_interfaces<I> ? coeff_deriv_patch2 : coeff_deriv_patch1;
        const double coefs[] = {-coeff_left, 1, -coeff_right};

        for (int dofs : {-1, 0, 1}) {
            if (0 <= I + dofs && I + dofs < n_inner_interfaces) {
                m_matrix.set_element(I, I + dofs, coefs[dofs + 1]);
            }
        }

        // Add coefficients for the periodic case.
        if constexpr (base_type::is_periodic && I == 0) {
            m_matrix.set_element(I, n_inner_interfaces - 1, coefs[0]);
        }
        if constexpr (base_type::is_periodic && I == n_inner_interfaces - 1) {
            m_matrix.set_element(I, 0, coefs[2]);
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

    /**
     * @brief Fill in the Ith line of the vector contening the orientation signs of the 
     * parallel grids. Keeping these signs helps for the computation of equivalent global
     * spline cross-derivatives.
     */
    template <std::size_t I>
    void set_sign_vector_line()
    {
        // If the perpendicular grids agree on one dimension (first or second).
        constexpr bool is_same_dimension
                = ((std::is_same_v<grid_perp1<I>, typename patch1<I>::Grid1>)
                   == (std::is_same_v<grid_perp2<I>, typename patch2<I>::Grid1>));

        // If the perpendicular grids agree on the direction.
        constexpr bool are_same_direction_perp = (extremity_1<I> != extremity_2<I>);

        // If we need to change the sign compared to the previous interface.
        constexpr bool change_sign = (are_same_direction_perp == is_same_dimension);
        constexpr int sign_changed = change_sign - !change_sign;

        m_vector_par_signs[I + 1] = m_vector_par_signs[I] * sign_changed;
    }

    /// @brief Set the vector C.
    template <typename IdxDerivs, class GridPar1D, std::size_t... I>
    void set_vector(
            IdxDerivs const& deriv_order,
            Idx<GridPar1D> const& slice_idx,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            std::integer_sequence<std::size_t, I...>)
    {
        // Check we start the loop with the first patch.
        static_assert(ddc::in_tags_v<
                      GridPar1D,
                      ddc::TypeSeq<
                              typename patch_first::Grid1,
                              typename patch_first::Grid2,
                              ddc::Deriv<typename patch_first::Dim1>,
                              ddc::Deriv<typename patch_first::Dim2>>>);

        // Locally transform the index into an integer to avoid type issues.
        IdxRange<GridPar1D> idx_range_1d_first_patch(m_idx_ranges.template get<patch_first>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (set_line_vector<I>(deriv_order, slice_idx_value, functions_and_derivs), ...);
    }

    /// @brief Fill in the Ith line of the vector C corresponding.
    template <std::size_t I, typename IdxDerivs>
    void set_line_vector(
            IdxDerivs const& deriv_order,
            int& slice_previous_idx_value,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        constexpr bool is_first_deriv = std::is_same_v<IdxDerivs, idx_deriv_perp>;

        Idx<grid_par1<I>> idx_to_slice_1;
        Idx<grid_par2<I>> idx_to_slice_2;
        base_type::template set_slice_indexes<
                I>(slice_previous_idx_value, idx_to_slice_1, idx_to_slice_2);

        // Compute the coefficient c_I (= approximated derivative) for the interface I.
        double lin_comb_funct = base_type::template get_approximated_deriv<
                I,
                is_first_deriv>(idx_to_slice_1, idx_to_slice_2, functions_and_derivs);

        /*
            If the orientation of the interface is not correct, we change the sign of the sum.
            Because c_{I 1|2} = - c_{I 2|1}.
        */
        int sign_sum = patch_order_matches_interfaces<I> - !patch_order_matches_interfaces<I>;
        if constexpr (!is_first_deriv) {
            /*
                > For the cross-derivative computation:
                In get_approx_cross_deriv(), the sign of the derivatives on Patch2 is changed
                if the edge orientation disagree. We need to correct the sign, if Patch1 was 
                ill-oriented. Because, 
                        d(f(-x)) = - df(-x).
            */
            sign_sum *= m_vector_par_signs
                    [ddc::type_seq_rank_v<grid_perp1<I>, typename base_type::Grid1DSeq>];
        }

        m_vector[I] = sign_sum * lin_comb_funct;

        // Add the boundary derivatives for global Hermite boundary conditions. ------------------
        constexpr bool is_lower_bound_deriv_dependent
                = (LowerBound == ddc::SplineBuilderClosure::HERMITE && I == 0);
        constexpr bool is_upper_bound_deriv_dependent
                = (UpperBound == ddc::SplineBuilderClosure::HERMITE && I == n_inner_interfaces - 1);

        if constexpr (is_lower_bound_deriv_dependent || is_upper_bound_deriv_dependent) {
            DerivFieldOnPatch_host<patch1<I>> function_and_derivs_1
                    = functions_and_derivs.template get<patch1<I>>();
            DerivFieldOnPatch_host<patch2<I>> function_and_derivs_2
                    = functions_and_derivs.template get<patch2<I>>();

            // Get the indices at the other bound of the interface.
            Idx<grid_perp1<I>> idx_perp_1 = this->template get_idx_interface<false, I, patch1<I>>();
            Idx<grid_perp2<I>> idx_perp_2 = this->template get_idx_interface<false, I, patch2<I>>();

            auto const [sign_1, sign_2] = get_sign_derivs<I, is_first_deriv>();

            typename base_type::template idx_deriv_perp1<I>
                    idx_deriv_bound_1(Idx<deriv_perp1<I>>(1), idx_perp_1, idx_to_slice_1);
            typename base_type::template idx_deriv_perp2<I>
                    idx_deriv_bound_2(Idx<deriv_perp2<I>>(1), idx_perp_2, idx_to_slice_2);

            double deriv_bound_1;
            double deriv_bound_2;
            if constexpr (is_first_deriv) {
                // Select the first derivatives at the other interfaces/bounds.
                deriv_bound_1 = function_and_derivs_1(idx_deriv_bound_1) * sign_1;
                deriv_bound_2 = function_and_derivs_2(idx_deriv_bound_2) * sign_2;
            } else {
                // Select the cross-derivatives at the other interfaces/bounds.
                typename base_type::template idx_cross_deriv1<I>
                        idx_cross_deriv_1(Idx<deriv_par1<I>>(1), idx_deriv_bound_1);
                typename base_type::template idx_cross_deriv2<I>
                        idx_cross_deriv_2(Idx<deriv_par2<I>>(1), idx_deriv_bound_2);

                deriv_bound_1 = function_and_derivs_1(idx_cross_deriv_1) * sign_1;
                deriv_bound_2 = function_and_derivs_2(idx_cross_deriv_2) * sign_2;
            }

            const double coeff_deriv_1 = m_deriv_coeffs_collect.template get<stored_interface<I>>()
                                                 .get_coeff_deriv_patch_1();
            const double coeff_deriv_2 = m_deriv_coeffs_collect.template get<stored_interface<I>>()
                                                 .get_coeff_deriv_patch_2();

            if constexpr (is_lower_bound_deriv_dependent) {
                m_vector[I] += patch_order_matches_interfaces<I> ? coeff_deriv_1 * deriv_bound_1
                                                                 : coeff_deriv_2 * deriv_bound_2;
            } else {
                m_vector[I] += patch_order_matches_interfaces<I> ? coeff_deriv_2 * deriv_bound_2
                                                                 : coeff_deriv_1 * deriv_bound_1;
            }
        }
    }

    /// @brief Once the interface derivatives computed, fill in the correct derivative fields.
    template <typename IdxDerivs, class GridPar1D, std::size_t... I>
    void update_derivatives(
            IdxDerivs const& deriv_order,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            Idx<GridPar1D> const& slice_idx,
            std::integer_sequence<std::size_t, I...>)
    {
        static_assert(ddc::in_tags_v<
                      GridPar1D,
                      ddc::TypeSeq<
                              typename patch_first::Grid1,
                              typename patch_first::Grid2,
                              ddc::Deriv<typename patch_first::Dim1>,
                              ddc::Deriv<typename patch_first::Dim2>>>);

        // Locally transform the index into an integer to avoid type issues.
        IdxRange<GridPar1D> idx_range_1d_first_patch(m_idx_ranges.template get<patch_first>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (update_derivatives_at_interface<I>(deriv_order, functions_and_derivs, slice_idx_value),
         ...);
    }

    /// @brief Associate the Ith derivative values to the correct first derivative or cross-derivatives.
    template <std::size_t I, typename IdxDerivs>
    void update_derivatives_at_interface(
            IdxDerivs const& deriv_order,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& slice_previous_idx_value)
    {
        constexpr bool is_first_deriv = std::is_same_v<IdxDerivs, idx_deriv_perp>;

        Idx<grid_par1<I>> idx_to_slice_1;
        Idx<grid_par2<I>> idx_to_slice_2;
        base_type::template set_slice_indexes<
                I>(slice_previous_idx_value, idx_to_slice_1, idx_to_slice_2);

        auto const [sign_1, sign_2] = get_sign_derivs<I, is_first_deriv>();

        // Get the indices at the bound of the interface.
        Idx<grid_perp1<I>> idx_perp_1 = this->template get_idx_interface<true, I, patch1<I>>();
        Idx<grid_perp2<I>> idx_perp_2 = this->template get_idx_interface<true, I, patch2<I>>();

        typename base_type::template idx_deriv_perp1<I>
                idx_interface_deriv_1(Idx<deriv_perp1<I>>(1), idx_perp_1, idx_to_slice_1);
        typename base_type::template idx_deriv_perp2<I>
                idx_interface_deriv_2(Idx<deriv_perp2<I>>(1), idx_perp_2, idx_to_slice_2);
        if constexpr (is_first_deriv) {
            // Update the first derivatives.
            functions_and_derivs.template get<patch1<I>>()(idx_interface_deriv_1)
                    = m_interface_derivatives[I] * sign_1;
            functions_and_derivs.template get<patch2<I>>()(idx_interface_deriv_2)
                    = m_interface_derivatives[I] * sign_2;
        } else {
            // The slice indices has to be a point at a corner, i.e. index range boundaries.
            IdxRange<grid_par1<I>> idx_range_par_1(m_idx_ranges.template get<patch1<I>>());
            IdxRange<grid_par2<I>> idx_range_par_2(m_idx_ranges.template get<patch2<I>>());
            assert((idx_to_slice_1 == idx_range_par_1.front())
                   || (idx_to_slice_1 == idx_range_par_1.back()));
            assert((idx_to_slice_2 == idx_range_par_2.front())
                   || (idx_to_slice_2 == idx_range_par_2.back()));

            typename base_type::template idx_cross_deriv1<I>
                    idx_cross_deriv_1(Idx<deriv_par1<I>>(1), idx_interface_deriv_1);
            typename base_type::template idx_cross_deriv2<I>
                    idx_cross_deriv_2(Idx<deriv_par2<I>>(1), idx_interface_deriv_2);

            // Update the cross-derivative.
            functions_and_derivs.template get<patch1<I>>()(idx_cross_deriv_1)
                    = m_interface_derivatives[I] * sign_1;
            functions_and_derivs.template get<patch2<I>>()(idx_cross_deriv_2)
                    = m_interface_derivatives[I] * sign_2;
        }
    }

    /**
     * @brief Get the sign of the interface derivatives for two patches of the given interface.
     * Due to different orientations of the patch axes, the signs can be different from an 
     * equivalent global spline. 
     */
    template <std::size_t I, bool is_first_deriv, std::enable_if_t<is_first_deriv, bool> = true>
    std::pair<int, int> get_sign_derivs() const
    {
        // If the directions of the two perpendicular grids disagree, the sign is changed.
        constexpr bool are_same_direction_perp = (extremity_1<I> != extremity_2<I>);
        constexpr int sign_deriv_perp_2 = to_sign(are_same_direction_perp);

        // If Patch1 follows the orientation of the sorted 1D grid sequence.
        constexpr bool is_per_1_well_oriented
                = (patch_order_matches_interfaces<I> && (extremity_1<I> == Extremity::BACK))
                  || (!patch_order_matches_interfaces<I> && (extremity_1<I> == Extremity::FRONT));

        // Change the sign of the derivatives if Patch1 is ill-oriented.
        constexpr int sign_1 = to_sign(is_per_1_well_oriented);
        // Change in addition the sign of the derivative on Patch2, if the grid directions disagree.
        constexpr int sign_2 = sign_1 * sign_deriv_perp_2;

        return std::pair(sign_1, sign_2);
    };

    /**
     * @brief Get the sign of the interface cross-derivatives for two patches of the given interface.
     * Due to different orientations of the patch axes, the signs can be different from an 
     * equivalent global spline. 
     */
    template <std::size_t I, bool is_first_deriv, std::enable_if_t<!is_first_deriv, bool> = true>
    std::pair<int, int> get_sign_derivs() const
    {
        // If Patch1 follows the orientation of the sorted 1D grid sequence.
        constexpr bool is_per_1_well_oriented
                = (patch_order_matches_interfaces<I> && (extremity_1<I> == Extremity::BACK))
                  || (!patch_order_matches_interfaces<I> && (extremity_1<I> == Extremity::FRONT));
        constexpr int sign_deriv_perp_1 = is_per_1_well_oriented - !is_per_1_well_oriented;
        // Sign for the parallel axis of Patch1.
        const int sign_deriv_par_1 = m_vector_par_signs
                [ddc::type_seq_rank_v<grid_perp1<I>, typename base_type::Grid1DSeq>];

        // If the perpendicular grids of the two patches agree.
        constexpr bool are_same_direction_perp = (extremity_1<I> != extremity_2<I>);
        // If the parallel grids of the two patches agree.
        constexpr bool are_same_direction_par = stored_interface<I>::orientations_agree;

        constexpr int sign_deriv_perp_2 = to_sign(are_same_direction_perp);
        constexpr int sign_deriv_par_2 = to_sign(are_same_direction_par);

        /*
            For Patch1: 
            Change the sign of the cross-derivative if the perpendicular grid of Patch1 is ill-oriented.
            Change the sign of the cross-derivative if the parallel grid of Patch1 is ill-oriented. 
        */
        const int sign_1 = sign_deriv_perp_1 * sign_deriv_par_1;
        /*
            For Patch2:
            Give the same sign as for the cross-derivative of Patch1. 
            Change the sign of the cross-derivative if the orientation of the perpendicular grid 
            of Patch2 disagrees with Patch1.
            Change the sign of the cross-derivative if the orientation of the parallel grid 
            of Patch2 disagrees with Patch1.
        */
        const int sign_2 = sign_1 * sign_deriv_perp_2 * sign_deriv_par_2;

        return std::pair(sign_1, sign_2);
    };
};
