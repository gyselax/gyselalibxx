

# File constant\_extrapolation\_rules\_onion.hpp

[**File List**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**spline**](dir_729d943c83b6b5573a69e28a4db4673a.md) **>** [**constant\_extrapolation\_rules\_onion.hpp**](constant__extrapolation__rules__onion_8hpp.md)

[Go to the documentation of this file](constant__extrapolation__rules__onion_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "onion_patch_locator.hpp"
#include "utils_patch_locators.hpp"


template <class PatchLocator>
struct ConstantExtrapolationRuleOnion
{
private:
    static_assert(
            is_onion_patch_locator_v<PatchLocator>,
            "Extrapolation rule only defined for OnionPatchLocator.");

    using PatchOrdering = typename PatchLocator::PatchOrdering;

    static constexpr std::size_t n_patches = ddc::type_seq_size_v<PatchOrdering>;

    using MinRadiusPatch = ddc::type_seq_element_t<0, PatchOrdering>;
    using MaxRadiusPatch = ddc::type_seq_element_t<n_patches - 1, PatchOrdering>;

    using R_min = typename MinRadiusPatch::Dim1;
    using Theta_min = typename MinRadiusPatch::Dim2;

    using R_max = typename MaxRadiusPatch::Dim1;
    using Theta_max = typename MaxRadiusPatch::Dim2;

    ddc::ConstantExtrapolationRule<R_min, Theta_min> const bc_r_min;
    ddc::ConstantExtrapolationRule<R_max, Theta_max> const bc_r_max;

public:
    explicit ConstantExtrapolationRuleOnion(Coord<R_min> const& r_min, Coord<R_max> const& r_max)
        : bc_r_min(r_min)
        , bc_r_max(r_max)
    {
    }

    template <class... Dim, template <typename P> typename SplinesOnPatch, class... Patches>
    KOKKOS_FUNCTION double operator()(
            Coord<Dim...> const& coord_extrap,
            MultipatchField<SplinesOnPatch, Patches...> const& patches_splines,
            int const out_of_bounds_idx) const
    {
        assert((out_of_bounds_idx == PatchLocator::outside_rmin_domain)
               || (out_of_bounds_idx == PatchLocator::outside_rmax_domain));

        if (out_of_bounds_idx == PatchLocator::outside_rmin_domain) {
            SplinesOnPatch<MinRadiusPatch> const min_spline
                    = patches_splines.template get<MinRadiusPatch>();
            return bc_r_min(coord_extrap, min_spline);
        } else {
            SplinesOnPatch<MaxRadiusPatch> const max_spline
                    = patches_splines.template get<MaxRadiusPatch>();
            return bc_r_max(coord_extrap, max_spline);
        }
    }
};
```


