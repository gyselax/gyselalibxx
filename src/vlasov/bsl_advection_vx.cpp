#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/ChunkSpan>
#include <ddc/Coordinate>
#include <ddc/DiscreteDomain>
#include <ddc/discretization>

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

DSpanSpXVx BslAdvectionVx::operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
        const
{
    const IDomainX& x_dom = get_domain<IDimX>(allfdistribu);
    const IDomainVx& vx_dom = get_domain<IDimVx>(allfdistribu);
    const IDomainSp& sp_dom = get_domain<IDimSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DFieldVx feet_coords(vx_dom);

    InterpolatorVxProxy interpolator_vx = m_interpolator_vx.preallocate();

    // Compute efield = -dPhi/dx where Phi is the electric_potential
    Chunk<double, BSDomainX> elecpot_spline_coef(m_spline_x_builder.spline_domain());
    m_spline_x_builder(elecpot_spline_coef, electric_potential);
    DFieldX efield(x_dom);
    for (IndexX ix : x_dom) {
        efield(ix) = -m_spline_x_evaluator.deriv(to_real(ix), elecpot_spline_coef);
    }

    for (IndexSp isp : sp_dom) {
        double const sqrt_me_on_mspecies = std::sqrt(
                m_species_info.mass()(m_species_info.ielec()) / m_species_info.mass()(isp));

        for (IndexX ix : x_dom) {
            // compute the displacement
            double const dvx = m_species_info.charge()(isp) * sqrt_me_on_mspecies * dt * efield(ix);

            // compute the coordinates of the feet
            for (IndexVx iv : vx_dom) {
                feet_coords(iv) = to_real(iv) - dvx;
            }

            // build a spline representation of the data
            interpolator_vx(allfdistribu[isp][ix], feet_coords);
        }
    }

    return allfdistribu;
}
