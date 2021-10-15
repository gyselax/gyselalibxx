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

#include "splineadvectionvx.hpp"

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionVx::SplineAdvectionVx(
        SpeciesInformation const& species_info,
        SplineXBuilder const& spline_x_builder,
        SplineEvaluator<BSplinesX> const& spline_x_evaluator,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_spline_vx_builder(spline_vx_builder)
    , m_spline_vx_evaluator(spline_vx_evaluator)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
    , m_species_info(species_info)
{
}

DSpanSpXVx SplineAdvectionVx::operator()(
        DSpanSpXVx allfdistribu,
        DViewX electric_potential,
        double dt) const
{
    assert(get_domain<MeshVx>(allfdistribu) == m_spline_vx_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(allfdistribu);
    const MDomainVx& vx_dom = get_domain<MeshVx>(allfdistribu);
    const MDomainSp& sp_dom = get_domain<MeshSp>(allfdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DBlockVx feet_coords(vx_dom);

    // Construct a domain over the bounded basis and allocate memory on this support
    Block<double, BSDomainVx> spline_coef(m_spline_vx_builder.spline_domain());

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
            m_spline_vx_builder(
                    spline_coef,
                    allfdistribu[isp][ix],
                    &m_derivs_vxmin,
                    &m_derivs_vxmax);

            // evaluate the function at the feet using the spline
            m_spline_vx_evaluator(allfdistribu[isp][ix], feet_coords.cview(), spline_coef.cview());
        }
    }

    return allfdistribu;
}
