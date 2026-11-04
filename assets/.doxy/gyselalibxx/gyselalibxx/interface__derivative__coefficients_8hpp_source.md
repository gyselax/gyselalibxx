

# File interface\_derivative\_coefficients.hpp

[**File List**](files.md) **>** [**interface\_derivatives**](dir_d1bd52a3e76a422151eefdcc4e15c189.md) **>** [**interface\_derivative\_coefficients.hpp**](interface__derivative__coefficients_8hpp.md)

[Go to the documentation of this file](interface__derivative__coefficients_8hpp.md)


```C++
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
    using associated_interface = InterfaceType;

public:
    using base_type::base_type;

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

    template <class Layout1, class Layout2>
    inline double get_approx_cross_deriv(
            DConstField<IdxRange1DPerp_1, Kokkos::HostSpace, Layout1> const& derivs_1,
            DConstField<IdxRange1DPerp_2, Kokkos::HostSpace, Layout2> const& derivs_2) const
    {
        return base_type::get_approx_cross_deriv(
                derivs_1,
                derivs_2,
                (InterfaceType::orientations_agree) - !(InterfaceType::orientations_agree));
    }

    template <class Layout1, class Layout2>
    inline double get_approx_cross_deriv(
            DConstField<IdxRange1DPerp_2, Kokkos::HostSpace, Layout2> const& derivs_2,
            DConstField<IdxRange1DPerp_1, Kokkos::HostSpace, Layout1> const& derivs_1) const
            requires(!std::is_same_v<IdxRange1DPerp_2, IdxRange1DPerp_1>)
    {
        return base_type::get_approx_cross_deriv(
                derivs_1,
                derivs_2,
                (InterfaceType::orientations_agree) - !(InterfaceType::orientations_agree));
    }
};
```


