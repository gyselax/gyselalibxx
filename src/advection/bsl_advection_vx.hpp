// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <i_interpolator.hpp>
#include <species_info.hpp>

#include "iadvectionvx.hpp"

template <class Geometry, class DDimV>
class BslAdvectionVelocity : public IAdvectionVelocity<Geometry, DDimV>
{
    using DDimSp = typename Geometry::DDimSp;
    using DDimX = typename Geometry::template spatial_dim_for<DDimV>;
    using FdistribuDDom = typename Geometry::FdistribuDDom;
    using SpatialDDom = typename Geometry::SpatialDDom;
    static_assert(std::is_same_v<FdistribuDDom, ddc::DiscreteDomain<DDimSp, DDimX, DDimV>>);
    using DElemV = ddc::DiscreteElement<DDimV>;
    using DElemSp = ddc::DiscreteElement<DDimSp>;
    using CDimV = typename DDimV::continuous_dimension_type;

private:
    IPreallocatableInterpolator<DDimV> const& m_interpolator_v;

public:
    explicit BslAdvectionVelocity(IPreallocatableInterpolator<DDimV> const& interpolator_v)
        : m_interpolator_v(interpolator_v)
    {
    }

    ~BslAdvectionVelocity() override = default;

    ddc::ChunkSpan<double, FdistribuDDom> operator()(
            ddc::ChunkSpan<double, FdistribuDDom> const allfdistribu,
            ddc::ChunkSpan<const double, SpatialDDom> const electric_field,
            double const dt) const override
    {
        FdistribuDDom const dom = allfdistribu.domain();
        ddc::DiscreteDomain<DDimV> const v_dom = ddc::select<DDimV>(dom);
        ddc::DiscreteDomain<DDimSp> const sp_dom = ddc::select<DDimSp>(dom);

        // pre-allocate some memory to prevent allocation later in loop
        ddc::Chunk<ddc::Coordinate<CDimV>, ddc::DiscreteDomain<DDimV>> feet_coords(v_dom);
        ddc::Chunk<double, ddc::DiscreteDomain<DDimV>> contiguous_slice(v_dom);
        InterpolatorProxy<DDimV> const interpolator_v = m_interpolator_v.preallocate();

        SpatialDDom const spatial_dom(allfdistribu.domain());

        auto c_dom = ddc::remove_dims_of(
                ddc::remove_dims_of(ddc::remove_dims_of(dom, sp_dom), spatial_dom),
                v_dom);
        using DElemC = typename decltype(c_dom)::discrete_element_type;
        using DElemSpatial = typename SpatialDDom::discrete_element_type;

        ddc::for_each(c_dom, [&](DElemC const ic) {
            ddc::for_each(sp_dom, [&](DElemSp const isp) {
                double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));

                ddc::for_each(spatial_dom, [&](DElemSpatial const ix) {
                    // compute the displacement
                    double const dvx = charge(isp) * sqrt_me_on_mspecies * dt * electric_field(ix);

                    // compute the coordinates of the feet
                    ddc::for_each(v_dom, [&](DElemV const iv) {
                        feet_coords(iv) = ddc::Coordinate<CDimV>(ddc::coordinate(iv) - dvx);
                    });

                    // copy the slice in contiguous memory
                    ddc::deepcopy(contiguous_slice, allfdistribu[ic][isp][ix]);

                    // build a spline representation of the data
                    interpolator_v(contiguous_slice, feet_coords.span_cview());

                    // copy back
                    ddc::deepcopy(allfdistribu[ic][isp][ix], contiguous_slice);
                });
            });
        });

        return allfdistribu;
    }
};
