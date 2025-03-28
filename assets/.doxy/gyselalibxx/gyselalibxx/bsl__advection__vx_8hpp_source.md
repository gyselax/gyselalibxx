

# File bsl\_advection\_vx.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**bsl\_advection\_vx.hpp**](bsl__advection__vx_8hpp.md)

[Go to the documentation of this file](bsl__advection__vx_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "iadvectionvx.hpp"
#include "iinterpolator.hpp"
#include "species_info.hpp"

template <class Geometry, class GridV>
class BslAdvectionVelocity : public IAdvectionVelocity<Geometry, GridV>
{
    using IdxRangeFdistribu = typename Geometry::IdxRangeFdistribu;
    using IdxRangeSpatial = typename Geometry::IdxRangeSpatial;
    using IdxSpatial = typename IdxRangeSpatial::discrete_element_type;
    using IdxV = Idx<GridV>;
    using DimV = typename GridV::continuous_dimension_type;
    using IdxRangeSpaceVelocity
            = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species>;

private:
    using PreallocatableInterpolatorType = interpolator_on_idx_range_t<
            IPreallocatableInterpolator,
            GridV,
            IdxRangeSpaceVelocity>;
    using InterpolatorType
            = interpolator_on_idx_range_t<IInterpolator, GridV, IdxRangeSpaceVelocity>;
    PreallocatableInterpolatorType const& m_interpolator_v;

public:
    explicit BslAdvectionVelocity(PreallocatableInterpolatorType const& interpolator_v)
        : m_interpolator_v(interpolator_v)
    {
    }

    ~BslAdvectionVelocity() override = default;

    Field<double, IdxRangeFdistribu> operator()(
            Field<double, IdxRangeFdistribu> const allfdistribu,
            Field<const double, IdxRangeSpatial> const electric_field,
            double const dt) const override
    {
        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFdistribu, Species, GridV>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;


        Kokkos::Profiling::pushRegion("BslAdvectionVelocity");
        IdxRangeFdistribu const idx_range = get_idx_range(allfdistribu);
        IdxRange<GridV> const idx_range_v = ddc::select<GridV>(idx_range);
        IdxRange<Species> const idx_range_sp = ddc::select<Species>(idx_range);

        FieldMem<double, typename InterpolatorType::batched_derivs_idx_range_type> derivs_min(
                m_interpolator_v.batched_derivs_idx_range_xmin(
                        ddc::remove_dims_of<Species>(idx_range)));
        FieldMem<double, typename InterpolatorType::batched_derivs_idx_range_type> derivs_max(
                m_interpolator_v.batched_derivs_idx_range_xmax(
                        ddc::remove_dims_of<Species>(idx_range)));
        ddc::parallel_fill(derivs_min, 0.);
        ddc::parallel_fill(derivs_max, 0.);

        // pre-allocate some memory to prevent allocation later in loop
        IdxRangeSpaceVelocity batched_feet_idx_range(idx_range);
        FieldMem<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords_alloc(batched_feet_idx_range);
        Field<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords(get_field(feet_coords_alloc));
        std::unique_ptr<InterpolatorType> const interpolator_v_ptr = m_interpolator_v.preallocate();
        InterpolatorType const& interpolator_v = *interpolator_v_ptr;

        IdxRangeSpatial const idx_range_spatial(get_idx_range(allfdistribu));

        IdxRangeBatch batch_idx_range(idx_range);

        ddc::for_each(idx_range_sp, [&](IdxSp const isp) {
            double const charge_proxy
                    = charge(isp); // TODO: consider proper way to access charge from device
            double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(IdxBatch const ib) {
                        IdxSpatial const ix(ib);
                        // compute the displacement
                        double const dvx
                                = charge_proxy * sqrt_me_on_mspecies * dt * electric_field(ix);

                        // compute the coordinates of the feet
                        for (IdxV const iv : idx_range_v) {
                            feet_coords(iv, ib) = Coord<DimV>(ddc::coordinate(iv) - dvx);
                        }
                    });
            interpolator_v(
                    allfdistribu[isp],
                    get_const_field(feet_coords),
                    get_const_field(derivs_min),
                    get_const_field(derivs_max));
        });

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
```


