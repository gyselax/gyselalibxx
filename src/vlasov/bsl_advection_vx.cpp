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

BslAdvectionVx::BslAdvectionVx(
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        IPreallocatableInterpolatorVx const& interpolator_vx)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_interpolator_vx(interpolator_vx)
{
}

DSpanSpXVx BslAdvectionVx::operator()(
        DSpanSpXVx const allfdistribu,
        DViewX const electric_field,
        double const dt) const
{
    IDomainX const& x_dom = get_domain<IDimX>(allfdistribu);
    IDomainVx const& vx_dom = get_domain<IDimVx>(allfdistribu);
    IDomainSp const& sp_dom = get_domain<IDimSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DFieldVx feet_coords(vx_dom);

    InterpolatorVxProxy const interpolator_vx = m_interpolator_vx.preallocate();

    for_each(sp_dom, [&](IndexSp const isp) {
        double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));

        for_each(x_dom, [&](IndexX const ix) {
            // compute the displacement
            double const dvx = charge(isp) * sqrt_me_on_mspecies * dt * electric_field(ix);

            // compute the coordinates of the feet
            for_each(vx_dom, [&](IndexVx const iv) { feet_coords(iv) = coordinate(iv) - dvx; });

            // build a spline representation of the data
            interpolator_vx(allfdistribu[isp][ix], feet_coords);
        });
    });

    return allfdistribu;
}
