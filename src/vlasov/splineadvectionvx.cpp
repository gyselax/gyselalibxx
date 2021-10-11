#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

#include <sll/null_boundary_value.h>
#include <sll/spline_evaluator.h>

#include "splineadvectionvx.h"

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionVx::SplineAdvectionVx(const BSplinesVx& bspl, const SplineVxBuilder& spl_interp)
    : SplineAdvectionVx(bspl, spl_interp, NullBoundaryValue::value, NullBoundaryValue::value)
{
}

SplineAdvectionVx::SplineAdvectionVx(
        const BSplinesVx& bspl,
        const SplineVxBuilder& spl_interp,
        const BoundaryValue& bc_left,
        const BoundaryValue& bc_right)
    : m_vx_spline_basis(bspl)
    , m_spline_vx_builder(spl_interp)
    , m_spline_vx_evaluator(bspl, bc_left, bc_right)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
{
}

DSpanXVx SplineAdvectionVx::operator()(
        DSpanXVx fdistribu,
        DViewX efield,
        int charge_species,
        double sqrt_me_on_mspecies,
        double dt) const
{
    assert(get_domain<MeshVx>(fdistribu) == m_spline_vx_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(fdistribu);
    const MDomainVx& vx_dom = get_domain<MeshVx>(fdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DBlockVx feet_coords(vx_dom);

    // Construct a domain over the bounded basis and allocate memory on this support
    Block<double, BSDomainVx> spline_coef(m_spline_vx_builder.spline_domain());

    for (MCoordX ix : x_dom) {
        // compute the displacement
        const double dvx = charge_species * sqrt_me_on_mspecies * dt * efield(ix);

        // compute the coordinates of the feet
        for (MCoordVx iv : vx_dom) {
            feet_coords(iv) = vx_dom.to_real(iv) - dvx;
        }

        // build a spline representation of the data
        m_spline_vx_builder(spline_coef, fdistribu[ix], &m_derivs_vxmin, &m_derivs_vxmax);

        // evaluate the function at the feet using the spline
        m_spline_vx_evaluator(fdistribu[ix], feet_coords.cview(), spline_coef.cview());
    }

    return fdistribu;
}
