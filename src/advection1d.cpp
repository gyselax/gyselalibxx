#include <cassert>
#include <cmath>
#include "advection1d.hpp"
#include "selalib_CXX/splines/spline_1d.h"

Advection1D::Advection1D(const BSplines& bspl, const Spline_interpolator_1D& spl_interp,
            double dt,
            double (*bc_left)(double),
            double (*bc_right)(double))
    : m_bspl(bspl), m_spline_interpolator(spl_interp),
      m_bc_left(bc_left), m_bc_right(bc_right), dt(dt)
{
    if ( !m_bspl.periodic )
    {
        assert(bc_left != nullptr);
        assert(bc_right != nullptr);
    }
}

void Advection1D::step(View& current_values, View&& velocity) const
{
    View x = m_spline_interpolator.get_interp_points();
    assert(current_values.extent(0) == x.extent(0));
    assert(current_values.extent(0) == velocity.extent(0));
    std::vector<double> new_points_ptr(x.extent(0));
    View new_points(new_points_ptr.data(), x.extent(0));

    Spline_1D spline(&m_bspl);
    m_spline_interpolator.compute_interpolant(spline, current_values);

    double xmin = m_bspl.xmin;
    double xmax = m_bspl.xmax;

    if (m_bspl.periodic)
    {
        for (size_t i(0); i < current_values.extent(0); ++i)
        {
            new_points[i] = x[i] - velocity[i] * dt;
            if ( new_points[i] < xmin || new_points[i] > xmax)
            {
                new_points[i] -= std::floor( (new_points[i] - xmin)/(xmax-xmin) )*(xmax-xmin);
            }
        }
        spline.eval_array(new_points, current_values);
    }
    else
    {
        for (size_t i(0); i < current_values.extent(0); ++i)
        {
            double new_point = x[i] - velocity[i] * dt;
            if ( new_point < xmin)
            {
                current_values[i] = m_bc_left(new_point);
            }
            else if ( new_point > xmax)
            {
                current_values[i] = m_bc_right(new_point);
            }
            else
            {
                current_values[i] = spline.eval(new_point);
            }
        }
    }
}
