#include "spline_2d.h"
#include "bsplines_non_uniform.h"
#include "bsplines_uniform.h"
#include <cassert>
#include <iostream>


Spline_2D::Spline_2D(const BSplines& bspl1, const BSplines& bspl2)
    : bcoef_ptr(new double[(bspl1.degree + bspl1.ncells)*(bspl2.degree + bspl2.ncells)]),
      bcoef(bcoef_ptr.get(), bspl1.degree + bspl1.ncells, bspl2.degree + bspl2.ncells),
      bspl1(bspl1),
      bspl2(bspl2)
{}

bool Spline_2D::belongs_to_space(const BSplines& bspline1, const BSplines& bspline2) const
{
    return (&bspl1 == &bspline1 && &bspl2 == &bspline2);
}

template <class T1, class T2, bool deriv1, bool deriv2,
          typename std::enable_if<std::is_base_of<BSplines, T1>::value>::type* = nullptr,
          typename std::enable_if<std::is_base_of<BSplines, T2>::value>::type* = nullptr>
double Spline_2D::eval_intern(double x1, double x2, const T1& bspl1, const T2& bspl2, mdspan_1d& vals1, mdspan_1d& vals2) const
{
    int jmin1, jmin2;

    if constexpr (deriv1)
    {
        bspl1.eval_deriv(x1, vals1, jmin1);
    }
    else
    {
        bspl1.eval_basis(x1, vals1, jmin1);
    }
    if constexpr (deriv2)
    {
        bspl2.eval_deriv(x2, vals2, jmin2);
    }
    else
    {
        bspl2.eval_basis(x2, vals2, jmin2);
    }

    double y = 0.0;
    for (int i(0); i<bspl1.degree+1; ++i)
    {
        for (int j(0); j<bspl2.degree+1; ++j)
        {
            y += bcoef(jmin1 + i, jmin2 + j) * vals1(i) * vals2(j);
        }
    }
    return y;
}

double Spline_2D::eval(const double x1, const double x2) const
{
    return eval_deriv<false, false>(x1, x2);
}

template<bool deriv1, bool deriv2>
double Spline_2D::eval_deriv(const double x1, const double x2) const
{
    double values1[bspl1.degree+1];
    double values2[bspl2.degree+1];
    mdspan_1d vals1(values1, bspl1.degree+1);
    mdspan_1d vals2(values2, bspl2.degree+1);

    if (bspl1.uniform)
    {
        if (bspl2.uniform)
        {
            return eval_intern<BSplines_uniform,BSplines_uniform,deriv1,deriv2>(x1, x2,
                    static_cast<const BSplines_uniform&>(bspl1),
                    static_cast<const BSplines_uniform&>(bspl2),
                    vals1, vals2);
        }
        else
        {
            return eval_intern<BSplines_uniform, BSplines_non_uniform,deriv1,deriv2>(x1, x2,
                    static_cast<const BSplines_uniform&>(bspl1),
                    static_cast<const BSplines_non_uniform&>(bspl2),
                    vals1, vals2);
        }
    }
    else
    {
        if (bspl2.uniform)
        {
            return eval_intern<BSplines_non_uniform,BSplines_uniform,deriv1,deriv2>(x1, x2,
                    static_cast<const BSplines_non_uniform&>(bspl1),
                    static_cast<const BSplines_uniform&>(bspl2),
                    vals1, vals2);
        }
        else
        {
            return eval_intern<BSplines_non_uniform,BSplines_non_uniform,deriv1,deriv2>(x1, x2,
                    static_cast<const BSplines_non_uniform&>(bspl1),
                    static_cast<const BSplines_non_uniform&>(bspl2),
                    vals1, vals2);
        }
    }
}

void Spline_2D::eval_array(mdspan_2d const& x1, mdspan_2d const& x2, mdspan_2d& y) const
{
    assert( x1.extent(0) == y.extent(0) );
    assert( x1.extent(1) == y.extent(1) );
    assert( x2.extent(0) == y.extent(0) );
    assert( x2.extent(1) == y.extent(1) );

    for (int i(0); i<x1.extent(0); ++i)
    {
        for (int j(0); j<x1.extent(1); ++j)
        {
            y(i,j) = eval(x1(i,j), x2(i,j));
        }
    }
}

template <class T1, typename std::enable_if<std::is_base_of<BSplines, T1>::value>::type* = nullptr,
          class T2, typename std::enable_if<std::is_base_of<BSplines, T2>::value>::type* = nullptr>
void Spline_2D::eval_array_loop(mdspan_2d const& x1, mdspan_2d const& x2, mdspan_2d& y) const
{
    double values1[bspl1.degree+1];
    double values2[bspl2.degree+1];
    mdspan_1d vals1(values1, bspl1.degree+1);
    mdspan_1d vals2(values2, bspl2.degree+1);

    int jmin1, jmin2;

    const T1& l_bspl1 = static_cast<const T1&>(bspl1);
    const T2& l_bspl2 = static_cast<const T1&>(bspl2);

    for (int i(0); i<x1.extent(0); ++i)
    {
        for (int j(0); j<x1.extent(1); ++j)
        {
            y(i,j) = eval_intern<T1, T2, false, false>(x1(i,j), x2(i,j), l_bspl1, l_bspl2, vals1, vals2);
        }
    }
}

template<bool deriv1, bool deriv2>
void Spline_2D::eval_array_deriv(mdspan_2d const& x1, mdspan_2d const& x2, mdspan_2d& y) const
{
    assert( x1.extent(0) == y.extent(0) );
    assert( x1.extent(1) == y.extent(1) );
    assert( x2.extent(0) == y.extent(0) );
    assert( x2.extent(1) == y.extent(1) );

    for (int i(0); i<x1.extent(0); ++i)
    {
        for (int j(0); j<x1.extent(1); ++j)
        {
            y(i,j) = eval_deriv<deriv1,deriv2>(x1(i,j), x2(i,j));
        }
    }
}

void Spline_2D::integrate_dim(mdspan_1d& y, const int dim) const
{
    assert( dim >=0 and dim < 2 );
    assert( y.extent(0) == bcoef.extent(1 - dim) );

    const BSplines& bspline ((dim == 0) ? bspl1 : bspl2);

    double values[bcoef.extent(dim)];
    mdspan_1d vals(values, bcoef.extent(dim));

    bspline.integrals(vals);

    if (dim==0)
    {
        for (int i(0); i<y.extent(0); ++i)
        {
            y(i) = 0;
            for (int j(0); j<bcoef.extent(0); ++j)
            {
                y(i) += bcoef(j,i) * vals(j);
            }
        }
    }
    else
    {
        for (int i(0); i<y.extent(0); ++i)
        {
            y(i) = 0;
            for (int j(0); j<bcoef.extent(1); ++j)
            {
                y(i) += bcoef(i,j) * vals(j);
            }
        }
    }
}

double Spline_2D::integrate() const
{
    double int_values[bcoef.extent(0)];
    mdspan_1d int_vals(int_values, bcoef.extent(0));

    integrate_dim(int_vals, 1);

    double values[bcoef.extent(0)];
    mdspan_1d vals(values, bcoef.extent(0));

    bspl1.integrals(vals);

    double y = 0.0;
    for (int i(0); i<bcoef.extent(0); ++i)
    {
        y += int_vals(i) * vals(i);
    }
    return y;
}