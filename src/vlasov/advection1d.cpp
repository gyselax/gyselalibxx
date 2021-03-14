#include <cassert>
#include <cmath>

#include "advection1d.h"
#include "spline_1d.h"

using namespace std;
using namespace std::experimental;

Advection1D::Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp)
    : m_bspl(bspl)
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
    : m_bspl(bspl)
    , m_spline_interpolator(spl_interp)
    , m_bc_left(bc_left)
    , m_bc_right(bc_right)
{
}

void Advection1D::operator()(DBlockView2D& fdistribu, double dt) const
{
    View1D<double> x = m_spline_interpolator.get_interp_points();
    unique_ptr<double[]> new_points_ptr = make_unique<double[]>(x.extent(0));
    View1D<double> new_points(new_points_ptr.get(), x.extent(0));

    unique_ptr<double[]> current_values_ptr = make_unique<double[]>(x.extent(0));
    View1D<double> current_values(current_values_ptr.get(), x.extent(0));

    Spline_1D spline(m_bspl, m_bc_left, m_bc_right);
    for (size_t vi = 0; vi < fdistribu.extent(1); ++vi) {
        double velocity = fdistribu.domain().mesh()[1]({vi})[0];

        auto current_values_nc = subspan(fdistribu.raw_view(), all, vi);
        for (size_t ii = 0; ii < current_values_nc.extent(0); ++ii) {
            current_values[ii] = current_values_nc[ii];
        }

        m_spline_interpolator.compute_interpolant(spline, current_values);
        for (size_t i = 0; i < current_values.extent(0); ++i) {
            new_points[i] = x[i] - velocity * dt;
        }
        spline.eval_array(new_points, current_values);
    }
}
