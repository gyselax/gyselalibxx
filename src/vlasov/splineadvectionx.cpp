#include <cassert>
#include <iosfwd>

#include <block_spline.h>
#include <deepcopy.h>
#include <spline_evaluator.h>

#include "blockview.h"
#include "mdomain.h"
#include "null_boundary_value.h"
#include "product_mdomain.h"
#include "product_mesh.h"
#include "rcoord.h"
#include "splineadvectionx.h"
#include "taggedvector.h"

#include <experimental/mdspan>

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionX::SplineAdvectionX(const BSplinesX& bspl, const SplineXBuilder& spl_interp)
    : m_x_spline_basis(bspl)
    , m_spline_x_builder(spl_interp)
    , m_bc_left(NullBoundaryValue::value)
    , m_bc_right(NullBoundaryValue::value)
{
    assert(bspl.is_periodic());
}

SplineAdvectionX::SplineAdvectionX(
        const BSplinesX& bspl,
        const SplineXBuilder& spl_interp,
        const BoundaryValue& bc_left,
        const BoundaryValue& bc_right)
    : m_x_spline_basis(bspl)
    , m_spline_x_builder(spl_interp)
    , m_bc_left(bc_left)
    , m_bc_right(bc_right)
{
}

DBlockSpanXVx SplineAdvectionX::operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt)
        const
{
    assert(get_domain<MeshX>(fdistribu) == m_spline_x_builder.interpolation_domain());

    const MDomainX& x_dom = get_domain<MeshX>(fdistribu);
    const MDomainVx& v_dom = get_domain<MeshVx>(fdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    DBlockX feet_coords(x_dom);
    //BlockX<RCoordX> feet_coords(x_dom);
    DBlockX contiguous_slice(x_dom);

    Block<BSplinesX, double> spline(m_x_spline_basis);
    SplineEvaluator spline_evaluator(spline, m_bc_left, m_bc_right);

    for (MCoordVx vii : v_dom) {
        // compute the displacement
        const double dx = mass_ratio * dt * v_dom.mesh().to_real(vii);

        // compute the coordinates of the feet
        for (MCoordX xii : x_dom) {
            feet_coords(xii) = RCoordX(x_dom.mesh().to_real(xii) - dx);
        }

        // copy the slice in contiguous memory
        deepcopy(contiguous_slice, fdistribu[vii]);

        // build a spline representation of the data
        m_spline_x_builder(spline, contiguous_slice);

        // evaluate the function at the feet using the spline
        spline_evaluator(contiguous_slice);

        // copy back
        deepcopy(fdistribu[vii], contiguous_slice);
    }

    return fdistribu;
}
