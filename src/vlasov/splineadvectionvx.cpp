#include <cassert>
#include <iosfwd>

#include <experimental/mdspan>

#include <ddc/BlockSpan>
#include <ddc/MDomain>
#include <ddc/ProductMDomain>
#include <ddc/ProductMesh>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

#include <sll/block_spline.h>
#include <sll/null_boundary_value.h>
#include <sll/spline_evaluator.h>

#include "splineadvectionvx.h"

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionVx::SplineAdvectionVx(const BSplinesVx& bspl, const SplineVxBuilder& spl_interp)
    : m_vx_spline_basis(bspl)
    , m_spline_vx_builder(spl_interp)
    , m_bc_left(NullBoundaryValue::value)
    , m_bc_right(NullBoundaryValue::value)
{
    assert(bspl.is_periodic());
}

SplineAdvectionVx::SplineAdvectionVx(
        const BSplinesVx& bspl,
        const SplineVxBuilder& spl_interp,
        const BoundaryValue& bc_left,
        const BoundaryValue& bc_right)
    : m_vx_spline_basis(bspl)
    , m_spline_vx_builder(spl_interp)
    , m_bc_left(bc_left)
    , m_bc_right(bc_right)
{
}

DSpanXVx SplineAdvectionVx::operator()(
        DSpanXVx fdistribu,
        DViewX efield,
        double sqrt_me_on_mspecies,
        double dt) const
{
    // assert(get_domain<MeshVx>(fdistribu) == m_spline_v_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(fdistribu);
    const MDomainVx& vx_dom = get_domain<MeshVx>(fdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DBlockVx feet_coords(vx_dom);
    //BlockVx<RCoordVx> feet_coords(x_dom);

    Block<double, SplineVxBuilder::interpolation_domain_type> interpolated_fdistribu(
            m_spline_vx_builder.interpolation_domain());

    Block<double, BSplinesVx> spline(m_vx_spline_basis);
    SplineEvaluator spline_evaluator(spline, m_bc_left, m_bc_right);

    for (MCoordX xii : x_dom) {
        // compute the displacement
        const double dvx = sqrt_me_on_mspecies * dt * efield(xii);

        // compute the coordinates of the feet
        for (MCoordVx vii : vx_dom) {
            feet_coords(vii) = RCoordVx(vx_dom.to_real(vii) - dvx);
        }

        // some_interpolation(interpolated_fdistribu, fdistribu[xii]);

        // build a spline representation of the data
        m_spline_vx_builder(spline, interpolated_fdistribu);

        // evaluate the function at the feet using the spline
        spline_evaluator(fdistribu[xii], feet_coords.cview());
    }

    return fdistribu;
}
