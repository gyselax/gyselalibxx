#include <cassert>
#include <cmath>
#include <vector>

#include "advectionx.h"
#include "spline_1d.h"

using namespace std;
using namespace std::experimental;

AdvectionX::AdvectionX(const BSplines& bspl, const Spline_interpolator_1D& spl_interp)
    : m_x_spline_basis(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(NullBoundaryValue::value)
    , m_bc_right(NullBoundaryValue::value)
{
    assert(bspl.periodic);
}

AdvectionX::AdvectionX(
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

DBlockViewXVx& AdvectionX::operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt) const
{
    //TODO:     assert(get_domain<Dim::X>(fdistribu) == m_spline_interpolator.domain());

    const MDomainX& x_dom = get_domain<Dim::X>(fdistribu);
    const MDomainVx& v_dom = get_domain<Dim::Vx>(fdistribu);

    // pre-allocate some memory to prevent allocation later in loop
    //TODO:     BlockX<RCoordX> feet_coords(x_dom);
    DBlockX feet_coords(x_dom);
    DBlockX contiguous_slice(x_dom);
    //TODO:     SplineX spline(m_x_spline_basis);
    Spline_1D spline(m_x_spline_basis);

    for (MCoordVx vii : v_dom) {
        // compute the displacement
        const double dx = mass_ratio * dt * v_dom.to_real(vii);

        // compute the coordinates of the feet
        for (MCoordX xii : x_dom) {
            feet_coords(xii) = RCoordX(x_dom.to_real(xii) - dx);
        }

        // copy the slice in contiguous memory
        //TODO:         deepcopy(contiguous_slice, fdistribu.slice(vii));
        deepcopy(contiguous_slice, fdistribu.slice(all, vii.get<Dim::Vx>()));

        // build a spline representation of the data
        //TODO:         m_spline_interpolator(spline, contiguous_slice);
        m_spline_interpolator.compute_interpolant(spline, contiguous_slice.raw_view());

        // evaluate the function at the feet using the spline
        spline.eval_array(feet_coords.raw_view(), contiguous_slice.raw_view());
        //TODO:         spline.eval_at(contiguous_slice, feet_coords);

        // copy back
        //TODO:         deepcopy(fdistribu.slice(vii), contiguous_slice);
        deepcopy(fdistribu.slice(all, vii.get<Dim::Vx>()), contiguous_slice);
    }

    return fdistribu;
}
