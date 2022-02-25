// SPDX-License-Identifier: MIT

#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/ChunkSpan>
#include <ddc/Coordinate>
#include <ddc/DiscreteDomain>
#include <ddc/discretization>
#include <ddc/for_each>

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
        SpeciesInformation const& species_info,
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        IPreallocatableInterpolatorVx const& interpolator_vx)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_interpolator_vx(interpolator_vx)
    , m_species_info(species_info)
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
        double const sqrt_me_on_mspecies = std::sqrt(
                m_species_info.mass()(m_species_info.ielec()) / m_species_info.mass()(isp));

        for_each(x_dom, [&](IndexX const ix) {
            // compute the displacement
            double const dvx
                    = m_species_info.charge()(isp) * sqrt_me_on_mspecies * dt * electric_field(ix);

            // compute the coordinates of the feet
            for_each(vx_dom, [&](IndexVx const iv) { feet_coords(iv) = to_real(iv) - dvx; });

            // build a spline representation of the data
            interpolator_vx(allfdistribu[isp][ix], feet_coords);
        });
    });

    return allfdistribu;
}
