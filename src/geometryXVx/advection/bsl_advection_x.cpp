// SPDX-License-Identifier: MIT

#include <cassert>

#include <ddc/ddc.hpp>

#include <species_info.hpp>

#include "bsl_advection_x.hpp"
#include "i_interpolator_x.hpp"

using namespace std;

BslAdvectionX::BslAdvectionX(IPreallocatableInterpolatorX const& interpolator)
    : m_interpolator(interpolator)
{
}

DSpanSpXVx BslAdvectionX::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    IDomainX const& x_dom = ddc::get_domain<IDimX>(allfdistribu);
    IDomainVx const& v_dom = ddc::get_domain<IDimVx>(allfdistribu);
    IDomainSp const& sp_dom = ddc::get_domain<IDimSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    FieldX<CoordX> feet_coords(x_dom);
    DFieldX contiguous_slice(x_dom);
    InterpolatorXProxy const interpolator = m_interpolator.preallocate();

    ddc::for_each(sp_dom, [&](IndexSp const isp) {
        double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));

        ddc::for_each(v_dom, [&](IndexVx const iv) {
            // compute the displacement
            double const dx = sqrt_me_on_mspecies * dt * ddc::coordinate(iv);

            // compute the coordinates of the feet
            ddc::for_each(x_dom, [&](IndexX const ix) {
                feet_coords(ix) = ddc::Coordinate<RDimX>(ddc::coordinate(ix) - dx);
            });

            // copy the slice in contiguous memory
            ddc::deepcopy(contiguous_slice, allfdistribu[isp][iv]);

            // interpolate the function at the feet using the provided interpolator
            interpolator(contiguous_slice, feet_coords.span_cview());

            // copy back
            ddc::deepcopy(allfdistribu[isp][iv], contiguous_slice);
        });
    });

    return allfdistribu;
}