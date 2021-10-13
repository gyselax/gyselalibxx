#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

#include <sll/null_boundary_value.h>
#include <sll/spline_evaluator.h>

#include <species_info.hpp>

#include "splineadvectionvx.h"

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionVx::SplineAdvectionVx(
        SpeciesInformation const& species_info,
        const BSplinesX& bspl_x,
        const SplineXBuilder& spl_x_interp,
        const BSplinesVx& bspl_vx,
        const SplineVxBuilder& spl_vx_interp)
    : SplineAdvectionVx(
            species_info,
            bspl_x,
            spl_x_interp,
            bspl_vx,
            spl_vx_interp,
            NullBoundaryValue::value,
            NullBoundaryValue::value)
{
}

SplineAdvectionVx::SplineAdvectionVx(
        SpeciesInformation const& species_info,
        const BSplinesX& bspl_x,
        const SplineXBuilder& spl_x_interp,
        const BSplinesVx& bspl_vx,
        const SplineVxBuilder& spl_vx_interp,
        const BoundaryValue& bc_vx_left,
        const BoundaryValue& bc_vx_right)
    : m_spline_x_builder(spl_x_interp)
    , m_spline_x_evaluator(bspl_x, NullBoundaryValue::value, NullBoundaryValue::value)
    , m_spline_vx_builder(spl_vx_interp)
    , m_spline_vx_evaluator(bspl_vx, bc_vx_left, bc_vx_right)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
    , m_species_info(species_info)
{
}

DSpanSpXVx SplineAdvectionVx::operator()(DSpanSpXVx fdistribu, DViewX electric_potential, double dt)
        const
{
    assert(get_domain<MeshVx>(fdistribu) == m_spline_vx_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(fdistribu);
    const MDomainVx& vx_dom = get_domain<MeshVx>(fdistribu);
    const MDomainSp& sp_dom = get_domain<MeshSp>(fdistribu);

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
            m_spline_vx_builder(spline_coef, fdistribu[isp][ix], &m_derivs_vxmin, &m_derivs_vxmax);

            // evaluate the function at the feet using the spline
            m_spline_vx_evaluator(fdistribu[isp][ix], feet_coords.cview(), spline_coef.cview());
        }
    }

    return fdistribu;
}
