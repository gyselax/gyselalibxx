#include <cassert>
#include <iosfwd>

#include <deepcopy.h>

#include "deprecated/bsplines.h"
#include "deprecated/spline_1d.h"

#include "blockview.h"
#include "mcoord.h"
#include "mdomain.h"
#include "null_boundary_value.h"
#include "product_mdomain.h"
#include "product_mesh.h"
#include "rcoord.h"
#include "spline_builder_1d.h"
#include "splineadvectionx.h"
#include "taggedarray.h"

#include <experimental/mdspan>

class BoundaryValue;

using namespace std;
using namespace std::experimental;

SplineAdvectionX::SplineAdvectionX(
        const deprecated::BSplines& bspl,
        const deprecated::SplineBuilder1D& spl_interp)
    : m_x_spline_basis(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(NullBoundaryValue::value)
    , m_bc_right(NullBoundaryValue::value)
{
    assert(bspl.is_periodic());
}

SplineAdvectionX::SplineAdvectionX(
        const deprecated::BSplines& bspl,
        const deprecated::SplineBuilder1D& spl_interp,
        const BoundaryValue& bc_left,
        const BoundaryValue& bc_right)
    : m_x_spline_basis(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(bc_left)
    , m_bc_right(bc_right)
{
}

DBlockSpanXVx SplineAdvectionX::operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt)
        const
{
    //TODO: spline on mesh
    //assert(get_domain<Dim::X>(fdistribu) == m_spline_interpolator.domain());

    const MDomainX& x_dom = get_domain<MeshX>(fdistribu);
    const MDomainVx& v_dom = get_domain<MeshVx>(fdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    //TODO: spline on mesh
    DBlockX feet_coords(x_dom);
    //BlockX<RCoordX> feet_coords(x_dom);
    DBlockX contiguous_slice(x_dom);
    //TODO: spline on mesh
    deprecated::Spline1D spline(m_x_spline_basis);
    //SplineX spline(m_x_spline_basis);

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
        //TODO: spline on mesh
        m_spline_interpolator.compute_interpolant(spline, contiguous_slice.raw_view());
        //m_spline_interpolator(spline, contiguous_slice);

        // evaluate the function at the feet using the spline
        //TODO: spline on mesh
        spline.eval_array(feet_coords.raw_view(), contiguous_slice.raw_view());
        //spline.eval_at(contiguous_slice, feet_coords);

        // copy back
        deepcopy(fdistribu[vii], contiguous_slice);
    }

    return fdistribu;
}
