// SPDX-License-Identifier: MIT

#pragma once

#include "deriv_reconstruction_coeffs.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "types.hpp"

template <class InterfaceType>
class InterfaceDerivCoeffs
    : public DerivReconstructionCoeffs<
              typename InterfaceType::Edge1::perpendicular_grid,
              typename InterfaceType::Edge2::perpendicular_grid,
              find_grid_t<
                      typename InterfaceType::Edge1::perpendicular_grid::continuous_dimension_type,
                      ddc::detail::TypeSeq<
                              typename InterfaceType::Edge1::associated_patch::BSplines1,
                              typename InterfaceType::Edge1::associated_patch::BSplines2>>,
              find_grid_t<
                      typename InterfaceType::Edge2::perpendicular_grid::continuous_dimension_type,
                      ddc::detail::TypeSeq<
                              typename InterfaceType::Edge2::associated_patch::BSplines1,
                              typename InterfaceType::Edge2::associated_patch::BSplines2>>,
              InterfaceType::Edge1::extremity,
              InterfaceType::Edge2::extremity>
{
    using base_type = DerivReconstructionCoeffs<
            typename InterfaceType::Edge1::perpendicular_grid,
            typename InterfaceType::Edge2::perpendicular_grid,
            find_grid_t<
                    typename InterfaceType::Edge1::perpendicular_grid::continuous_dimension_type,
                    ddc::detail::TypeSeq<
                            typename InterfaceType::Edge1::associated_patch::BSplines1,
                            typename InterfaceType::Edge1::associated_patch::BSplines2>>,
            find_grid_t<
                    typename InterfaceType::Edge2::perpendicular_grid::continuous_dimension_type,
                    ddc::detail::TypeSeq<
                            typename InterfaceType::Edge2::associated_patch::BSplines1,
                            typename InterfaceType::Edge2::associated_patch::BSplines2>>,
            InterfaceType::Edge1::extremity,
            InterfaceType::Edge2::extremity>;

    using Patch1 = typename InterfaceType::Edge1::associated_patch;
    using Patch2 = typename InterfaceType::Edge2::associated_patch;

    using IdxRange1DPerp_1 = typename base_type::IdxRange1DPerp_1;
    using IdxRange1DPerp_2 = typename base_type::IdxRange1DPerp_2;

public:
    /// @brief Interface between the two involved patches.
    using associated_interface = InterfaceType;

public:
    using base_type::base_type;

    /**
     * @brief Get the coefficient (a) or (b) in front of the derivative on the given patch of the given Interface.
     * @return The value of the coefficient (a) or (b).
     */
    template <class Patch>
    double get_coeff_deriv_on_patch() const requires(is_patch_v<Patch>)
    {
        static_assert(
                (std::is_same_v<Patch, Patch1>) || (std::is_same_v<Patch, Patch2>),
                "The given Patch template parameter is not one of the Patch of the given "
                "Interface.");
        if constexpr (std::is_same_v<Patch, Patch1>) {
            return this->m_coeff_deriv_patch_1;
        } else {
            return this->m_coeff_deriv_patch_2;
        }
    }

    /**
     * @brief Get the linear combination of the derivatives values (c).
     *
     * @anchor get_approx_cross_deriv
     * This linear combination corresponds to an approximation of the
     * interface cross-derivative.
     *
     * It is more appropriate than get_approx_deriv to compute the
     * cross-derivatives. The computation is the same.
     * However, in addition, the sign of the derivatives on the patch 2 is
     * changed in the case where the orientations of the two patches do not
     * match.
     * E.g. dx(f(-x)) = - dxf(-x)
     *
     * @param derivs_1 Derivatives at the interpolation points on patch 1.
     * @param derivs_2 Derivatives at the interpolation points on patch 2.
     * @return the linear combination of the function values (c).
     */
    template <class Layout1, class Layout2>
    inline double get_approx_cross_deriv(
            DConstField<IdxRange1DPerp_1, Kokkos::HostSpace, Layout1> const& derivs_1,
            DConstField<IdxRange1DPerp_2, Kokkos::HostSpace, Layout2> const& derivs_2) const
    {
        int const sign_deriv_2
                = (InterfaceType::orientations_agree) - !(InterfaceType::orientations_agree);
        return this->get_approx_with_change_of_sign(derivs_1, derivs_2, sign_deriv_2);
    }

    /**
     * @brief Get the linear combination of the derivatives values (c).
     * See @ref get_approx_cross_deriv.
     * @param derivs_1 Derivatives at the interpolation points on patch 1.
     * @param derivs_2 Derivatives at the interpolation points on patch 2.
     * @return the linear combination of the function values (c).
     */
    template <class Layout1, class Layout2>
    inline double get_approx_cross_deriv(
            DConstField<IdxRange1DPerp_2, Kokkos::HostSpace, Layout2> const& derivs_2,
            DConstField<IdxRange1DPerp_1, Kokkos::HostSpace, Layout1> const& derivs_1) const
            requires(!std::is_same_v<IdxRange1DPerp_2, IdxRange1DPerp_1>)
    {
        int const sign_deriv_2
                = (InterfaceType::orientations_agree) - !(InterfaceType::orientations_agree);
        return this->get_approx_with_change_of_sign(derivs_1, derivs_2, sign_deriv_2);
    }
};
