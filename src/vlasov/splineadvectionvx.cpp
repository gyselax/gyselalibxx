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
        const BSplinesVx& bspl,
        const SplineVxBuilder& spl_interp)
    : SplineAdvectionVx(
            species_info,
            bspl,
            spl_interp,
            NullBoundaryValue::value,
            NullBoundaryValue::value)
{
}

SplineAdvectionVx::SplineAdvectionVx(
        SpeciesInformation const& species_info,
        const BSplinesVx& bspl,
        const SplineVxBuilder& spl_interp,
        const BoundaryValue& bc_left,
        const BoundaryValue& bc_right)
    : m_spline_vx_builder(spl_interp)
    , m_spline_vx_evaluator(bspl, bc_left, bc_right)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
    , m_species_info(species_info)
{
}

DSpanSpXVx SplineAdvectionVx::operator()(DSpanSpXVx fdistribu, DViewX efield, double dt) const
{
    assert(get_domain<MeshVx>(fdistribu) == m_spline_vx_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(fdistribu);
    const MDomainVx& vx_dom = get_domain<MeshVx>(fdistribu);
    const MDomainSp& sp_dom = get_domain<MeshSp>(fdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DBlockVx feet_coords(vx_dom);

    // Construct a domain over the bounded basis and allocate memory on this support
    Block<double, BSDomainVx> spline_coef(m_spline_vx_builder.spline_domain());

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
