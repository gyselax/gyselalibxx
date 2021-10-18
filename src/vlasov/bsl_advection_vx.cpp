#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

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
    const MDomainX& x_dom = get_domain<MeshX>(allfdistribu);
    const MDomainVx& vx_dom = get_domain<MeshVx>(allfdistribu);
    const MDomainSp& sp_dom = get_domain<MeshSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DBlockVx feet_coords(vx_dom);

    InterpolatorVxProxy interpolator_vx = m_interpolator_vx.preallocate();

    // Compute efield = -dPhi/dx where Phi is the electric_potential
    Block<double, BSDomainX> elecpot_spline_coef(m_spline_x_builder.spline_domain());
    m_spline_x_builder(elecpot_spline_coef, electric_potential);
    DBlockX efield(x_dom);
    for (MCoordX ix : x_dom) {
        efield(ix) = -m_spline_x_evaluator.deriv(x_dom.to_real(ix), elecpot_spline_coef.cview());
    }

    for (MCoordSp isp : sp_dom) {
        const double sqrt_me_on_mspecies = std::sqrt(
                m_species_info.mass()(m_species_info.ielec()) / m_species_info.mass()(isp));

        for (MCoordX ix : x_dom) {
            // compute the displacement
            const double dvx = m_species_info.charge()(isp) * sqrt_me_on_mspecies * dt * efield(ix);

            // compute the coordinates of the feet
            for (MCoordVx iv : vx_dom) {
                feet_coords(iv) = vx_dom.to_real(iv) - dvx;
            }

            // build a spline representation of the data
            interpolator_vx(allfdistribu[isp][ix], feet_coords.cview());
        }
    }

    return allfdistribu;
}
