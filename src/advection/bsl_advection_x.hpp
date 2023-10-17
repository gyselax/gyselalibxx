// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <i_interpolator.hpp>
#include <species_info.hpp>

#include "iadvectionx.hpp"

template <class Geometry, class DDimX>
class BslAdvectionSpatial : public IAdvectionSpatial<Geometry, DDimX>
{
    using DDimSp = typename Geometry::DDimSp;
    using DDimV = typename Geometry::template velocity_dim_for<DDimX>;
    using DDom = typename Geometry::FdistribuDDom;
    using DElemX = ddc::DiscreteElement<DDimX>;
    using DElemV = ddc::DiscreteElement<DDimV>;
    using DElemSp = ddc::DiscreteElement<DDimSp>;
    using DElemSpV = ddc::DiscreteElement<DDimSp, DDimV>;
    using CDimX = typename DDimX::continuous_dimension_type;

private:
    IPreallocatableInterpolator<DDimX> const& m_interpolator_x;

public:
    BslAdvectionSpatial(IPreallocatableInterpolator<DDimX> const& interpolator_x)
        : m_interpolator_x(interpolator_x)
    {
    }

    ~BslAdvectionSpatial() override = default;

    ddc::ChunkSpan<double, DDom> operator()(
            ddc::ChunkSpan<double, DDom> const allfdistribu,
            double const dt) const override
    {
        Kokkos::Profiling::pushRegion("BslAdvectionSpatial");
        DDom const dom = allfdistribu.domain();
        ddc::DiscreteDomain<DDimX> const x_dom = ddc::select<DDimX>(dom);
        ddc::DiscreteDomain<DDimV> const v_dom = ddc::select<DDimV>(dom);
        ddc::DiscreteDomain<DDimSp> const sp_dom = ddc::select<DDimSp>(dom);

        // pre-allocate some memory to prevent allocation later in loop
        ddc::Chunk<ddc::Coordinate<CDimX>, ddc::DiscreteDomain<DDimX>> feet_coords(x_dom);
        ddc::Chunk<double, ddc::DiscreteDomain<DDimX>> contiguous_slice(x_dom);
        std::unique_ptr<IInterpolator<DDimX>> const interpolator_x_ptr
                = m_interpolator_x.preallocate();
        IInterpolator<DDimX> const& interpolator_x = *interpolator_x_ptr;

        auto c_dom = ddc::remove_dims_of(
                dom,
                ddc::DiscreteDomain<DDimSp, DDimX, DDimV>(sp_dom, x_dom, v_dom));
        using DElemC = typename decltype(c_dom)::discrete_element_type;

        ddc::for_each(c_dom, [&](DElemC const ic) {
            ddc::for_each(sp_dom, [&](DElemSp const isp) {
                double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));

                ddc::for_each(v_dom, [&](DElemV const iv) {
                    // compute the displacement
                    double const dx = sqrt_me_on_mspecies * dt * ddc::coordinate(iv);

                    // compute the coordinates of the feet
                    ddc::for_each(x_dom, [&](DElemX const ix) {
                        feet_coords(ix) = ddc::Coordinate<CDimX>(ddc::coordinate(ix) - dx);
                    });

                    // copy the slice in contiguous memory
                    ddc::deepcopy(contiguous_slice, allfdistribu[ic][isp][iv]);

                    // interpolate the function at the feet using the provided interpolator
                    interpolator_x(contiguous_slice, feet_coords.span_cview());

                    // copy back
                    ddc::deepcopy(allfdistribu[ic][isp][iv], contiguous_slice);
                });
            });
        });

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
