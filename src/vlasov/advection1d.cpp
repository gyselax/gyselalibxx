#include <cassert>
#include <cmath>
#include <vector>

#include "advection1d.h"
#include "spline_1d.h"

using namespace std;
using namespace std::experimental;

Advection1D::Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp)
    : m_x_spline_basis(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(NullBoundaryValue::value)
    , m_bc_right(NullBoundaryValue::value)
{
    assert(bspl.periodic);
}

Advection1D::Advection1D(
        const BSplines& bspl,
        const Spline_interpolator_1D& spl_interp,
        const BoundaryValue& bc_left,
        const BoundaryValue& bc_right)
    : m_x_spline_basis(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(bc_left)
    , m_bc_right(bc_right)
{
}

DBlockViewXVx& Advection1D::operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt) const
{
    assert(get_domain<Dim::X>(fdistribu) == m_spline_interpolator.domain());

    const MDomainX& x_dom = get_domain<Dim::X>(fdistribu);
    const MDomainVx& v_dom = get_domain<Dim::Vx>(fdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    BlockX<RCoordX> feet_coords = allocate_bloc<RCoordX>(x_dom);
    DBlockX contiguous_slice = allocate_dbloc(x_dom);
    SplineX spline = allocate_spline(m_x_spline_basis);

    for (MCoordVx vii : v_dom) {
        // compute the displacement
        const double dx = mass_ratio * dt * v_dom.to_real(vii);

        // compute the coordinates of the feet
        for (MCoordX xii : x_dom) {
            feet_coords(xii) = x_dom.to_real(xii) - dx;
        }

        // copy the slice in contiguous memory
        deepcopy(contiguous_slice, fdistribu.slice(vii));

        // build a spline representation of the data
        m_spline_interpolator(spline, contiguous_slice);

        // evaluate the function at the feet using the spline
        spline.eval_at(contiguous_slice, feet_coords);

        // copy back
        deepcopy(fdistribu.slice(vii), contiguous_slice);
    }

    return fdistribu;
}
