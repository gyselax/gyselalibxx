// SPDX-License-Identifier: MIT

#pragma once

#include <typeinfo>

#include <ddc/ddc.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "derivative_calculator_common.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "interface_derivative_coefficients_collection.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"
#include "view.hpp"

template <class Connectivity, class Grid1D, class PatchSeq, class DerivCoeffCollection>
class InterfacesDerivativeCalculator;

/**
  * @brief Class to compute the interface derivatives along a given direction 
  * on all the interfaces with an approximation formula.
  * 
  * This operator is implemented for conforming equivalent global meshes. 
  * 
  * When we call the operator .solve_deriv() or .solve_cross_deriv(), the operator loops
  * over the interface to compute the vector C_trunc to determine all the interface derivatives 
  * and update the DerivField given in input. 
  * 
  * The vector C_trunc is given by the coefficients (see README),
  * @f$ \{\sum_{k = - N_{reduc}}^{N_{reduc}} \omega_{k, N_{reduc}, N_{reduc}}^{i_I} f_{i_I+k} \}_I @f$
  * 
  * with 
  *     * @f$ N_{reduc} @f$ the number of chosen cells for the approximation, 
  *     * @f$ i_I @f$ the index of the grid at the Ith interface, 
  *     * @f$ f_{i_I+k} @f$ the function values if we compute the first derivatives, or the first 
  * derivatives if we compute the cross-derivatives. 
  * 
  * @tparam Connectivity A MultipatchConnectivity class describing all the patch connections.
  * @tparam Grid1D A given direction.
  * @tparam Patches List of patches containing all the involved patches in the given direction. 
  * @tparam DerivCoeffCollection A InterfaceDerivCoeffsCollection
  * that stores the InterfaceDerivCoeffs needed to compute the interface derivatives. 
  */
template <class Connectivity, class Grid1D, class... Patches, class DerivCoeffCollection>
class InterfacesDerivativeCalculator<
        Connectivity,
        Grid1D,
        ddc::TypeSeq<Patches...>,
        DerivCoeffCollection>
    : public DerivativeCalculatorCommon<
              InterfacesDerivativeCalculator<
                      Connectivity,
                      Grid1D,
                      ddc::TypeSeq<Patches...>,
                      DerivCoeffCollection>,
              Connectivity,
              Grid1D,
              ddc::TypeSeq<Patches...>,
              DerivCoeffCollection>
{
    using base_type = DerivativeCalculatorCommon<
            InterfacesDerivativeCalculator<
                    Connectivity,
                    Grid1D,
                    ddc::TypeSeq<Patches...>,
                    DerivCoeffCollection>,
            Connectivity,
            Grid1D,
            ddc::TypeSeq<Patches...>,
            DerivCoeffCollection>;

    friend base_type;

    using typename base_type::eval_cross_deriv;
    using typename base_type::eval_deriv;

private:
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
    using deriv_par1 = typename base_type::template deriv_par1<I>;
    template <std::size_t I>
    using deriv_par2 = typename base_type::template deriv_par2<I>;

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using patch_first = typename base_type::patch_first;

public:
    /**
     * @brief Instantiate InterfacesDerivativeCalculator. 
     *  
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param deriv_coeffs_collect InterfaceDerivCoeffsCollection containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     */
    InterfacesDerivativeCalculator(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            DerivCoeffCollection const& deriv_coeffs_collect)
        : base_type(idx_ranges, deriv_coeffs_collect)
    {
    }

private:
    template <typename eval_type, class IdxPar>
    void solve(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        // Check the input types.
        static_assert(
                (std::is_same_v<IdxPar, typename patch_first::Idx1>)
                        || (std::is_same_v<IdxPar, typename patch_first::Idx2>),
                "The given index has to be a 1D index defined on the first patch. See "
                "InterfacesDerivativeCalculator<...>::patch_first.");

        // Update derivatives.
        update_derivatives<eval_type>(
                functions_and_derivs,
                idx_par,
                std::make_integer_sequence<std::size_t, base_type::n_inner_interfaces> {});
    }


    /// @brief Once the interface derivatives computed, fill in the correct derivative fields.
    template <typename eval_type, class GridPar1D, std::size_t... I>
    void update_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            Idx<GridPar1D> const& idx_to_slice,
            std::integer_sequence<std::size_t, I...>)
    {
        static_assert(
                (std::is_same_v<GridPar1D, typename patch_first::Grid1>)
                || (std::is_same_v<GridPar1D, typename patch_first::Grid2>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename patch_first::Dim1>>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename patch_first::Dim2>>));

        // Locally transform the index into an integer to avoid type issues.
        IdxRange<GridPar1D> idx_range_1d_first_patch(
                this->m_idx_ranges.template get<patch_first>());
        int value_idx_to_slice = (idx_to_slice - idx_range_1d_first_patch.front()).value();
        (update_derivatives_at_interface<I, eval_type>(functions_and_derivs, value_idx_to_slice),
         ...);
    }


    /// @brief Associate the Ith derivative values to the correct first derivative or cross-derivative.
    template <std::size_t I, typename eval_type>
    void update_derivatives_at_interface(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& value_previous_idx_to_slice)
    {
        Idx<grid_par1<I>> idx_to_slice_1;
        Idx<grid_par2<I>> idx_to_slice_2;
        base_type::template set_slice_indexes<
                I>(value_previous_idx_to_slice, idx_to_slice_1, idx_to_slice_2);

        double const reconstructed_deriv = base_type::template get_approximated_deriv<
                I,
                eval_type>(idx_to_slice_1, idx_to_slice_2, functions_and_derivs);

        // Get the correct indices for the slices.
        Idx<grid_perp1<I>> idx_perp_1 = this->template get_idx_interface<true, I, patch1<I>>();
        Idx<grid_perp2<I>> idx_perp_2 = this->template get_idx_interface<true, I, patch2<I>>();

        std::array<const int, 2> signs = get_signs<I, eval_type>();

        // Update the derivatives or cross-derivatives.
        Idx<typename base_type::template deriv_perp1<I>> i_dperp1(1);
        Idx<typename base_type::template deriv_perp2<I>> i_dperp2(1);
        if constexpr (std::is_same_v<eval_type, eval_deriv>) {
            typename base_type::template idx_deriv_perp1<I>
                    idx_deriv1(i_dperp1, idx_perp_1, idx_to_slice_1);
            typename base_type::template idx_deriv_perp2<I>
                    idx_deriv2(i_dperp2, idx_perp_2, idx_to_slice_2);

            functions_and_derivs.template get<patch1<I>>()(idx_deriv1)
                    = reconstructed_deriv * signs[0];
            functions_and_derivs.template get<patch2<I>>()(idx_deriv2)
                    = reconstructed_deriv * signs[1];
        } else {
            typename base_type::template idx_cross_deriv1<I>
                    idx_deriv1(i_dperp1, Idx<deriv_par1<I>>(1), idx_perp_1, idx_to_slice_1);
            typename base_type::template idx_cross_deriv2<I>
                    idx_deriv2(i_dperp2, Idx<deriv_par2<I>>(1), idx_perp_2, idx_to_slice_2);

            functions_and_derivs.template get<patch1<I>>()(idx_deriv1)
                    = reconstructed_deriv * signs[0];
            functions_and_derivs.template get<patch2<I>>()(idx_deriv2)
                    = reconstructed_deriv * signs[1];
        }
    }

    /**
     * @brief Get the signs of the correction on the derivatives. 
     * The sign of the computed interface derivative can be wrong in the case
     * the axes do not agree. 
     */
    template <std::size_t I, typename eval_type>
    static std::array<const int, 2> get_signs()
    {
        // If the parallel grids of the two patches agree.
        constexpr bool are_same_direction_par
                = std::is_same_v<eval_type, eval_deriv>
                          ? true
                          : base_type::template interface<I>::orientations_agree;

        constexpr bool change_sign_1 = (base_type::template extremity_1<I> == Extremity::FRONT);
        constexpr bool change_sign_2 = ((base_type::template extremity_2<I> == Extremity::BACK)
                                        && are_same_direction_par)
                                       || ((base_type::template extremity_2<I> == Extremity::FRONT)
                                           && !are_same_direction_par);

        constexpr int sign_1 = !change_sign_1 - change_sign_1;
        constexpr int sign_2 = !change_sign_2 - change_sign_2;
        return std::array<const int, 2> {sign_1, sign_2};
    }
};
