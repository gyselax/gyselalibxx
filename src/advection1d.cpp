#include "advection1d.hpp"

#include <cassert>
#include <cmath>

#include "selalib_CXX/splines/spline_1d.h"

Advection1D::Advection1D(const BSplines& bspl,
    const Spline_interpolator_1D& spl_interp)
    : m_bspl(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(NullBoundaryValue::value)
    , m_bc_right(NullBoundaryValue::value)
{
    assert(bspl.periodic);
}

Advection1D::Advection1D(const BSplines& bspl,
    const Spline_interpolator_1D& spl_interp,
    const BoundaryValue& bc_left,
    const BoundaryValue& bc_right)
    : m_bspl(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(bc_left)
    , m_bc_right(bc_right)
{
}

void Advection1D::operator()(View& current_values, View&& velocity, double dt) const
{
    View x = m_spline_interpolator.get_interp_points();
    assert(current_values.extent(0) == x.extent(0));
    assert(current_values.extent(0) == velocity.extent(0));
    std::vector<double> new_points_ptr(x.extent(0));
    View new_points(new_points_ptr.data(), x.extent(0));

    Spline_1D spline(m_bspl, m_bc_left, m_bc_right);
    m_spline_interpolator.compute_interpolant(spline, current_values);

    double xmin = m_bspl.xmin;
    double xmax = m_bspl.xmax;

    for (size_t i(0); i < current_values.extent(0); ++i) {
        new_points[i] = x[i] - velocity[i] * dt;
    }
    spline.eval_array(new_points, current_values);
}
