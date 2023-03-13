// SPDX-License-Identifier: MIT

#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/null_boundary_value.hpp>
#include <sll/spline_evaluator.hpp>

#include <species_info.hpp>

#include "bsl_advection_vx.hpp"
#include "i_interpolator_vx.hpp"
#include "i_interpolator_x.hpp"

class BoundaryValue;

using namespace std;
using namespace std::experimental;

BslAdvectionVx::BslAdvectionVx(IPreallocatableInterpolatorVx const& interpolator_vx)
    : m_interpolator_vx(interpolator_vx)
{
}

DSpanSpXVx BslAdvectionVx::operator()(
        DSpanSpXVx const allfdistribu,
        DViewX const electric_field,
        double const dt) const
{
    IDomainX const& x_dom = ddc::get_domain<IDimX>(allfdistribu);
    IDomainVx const& vx_dom = ddc::get_domain<IDimVx>(allfdistribu);
    IDomainSp const& sp_dom = ddc::get_domain<IDimSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    FieldVx<CoordVx> feet_coords(vx_dom);

    InterpolatorVxProxy const interpolator_vx = m_interpolator_vx.preallocate();

    ddc::for_each(sp_dom, [&](IndexSp const isp) {
        double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));

        ddc::for_each(x_dom, [&](IndexX const ix) {
            // compute the displacement
            double const dvx = charge(isp) * sqrt_me_on_mspecies * dt * electric_field(ix);

            // compute the coordinates of the feet
            ddc::for_each(vx_dom, [&](IndexVx const iv) {
                feet_coords(iv) = ddc::Coordinate<RDimVx>(ddc::coordinate(iv) - dvx);
            });

            // build a spline representation of the data
            interpolator_vx(allfdistribu[isp][ix], feet_coords);
        });
    });

    return allfdistribu;
}
