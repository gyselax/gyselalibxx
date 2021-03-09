#include <cassert>
#include <type_traits>
#include "spline_1d.h"
#include "bsplines_non_uniform.h"
#include "bsplines_uniform.h"


Spline_1D::Spline_1D(const BSplines& bspl)
    : bcoef_ptr(new double[bspl.degree + bspl.ncells]),
      bcoef(bcoef_ptr.get(), bspl.degree + bspl.ncells),
      bspl(bspl)
{}

bool Spline_1D::belongs_to_space(const BSplines& bspline) const
{
    return &bspl == &bspline;
}

template <class T, typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
double Spline_1D::eval_intern(double x, const T& bspl, mdspan_1d& vals) const
{
    int jmin;

    bspl.eval_basis(x, vals, jmin);

    double y = 0.0;
    for (int i(0); i<bspl.degree+1; ++i)
    {
        y += bcoef(jmin+i) * vals(i);
    }
    return y;
}

template <class T, typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
double Spline_1D::eval_deriv_intern(double x, const T& bspl, mdspan_1d& vals) const
{
    int jmin;

    bspl.eval_deriv(x, vals, jmin);

    double y = 0.0;
    for (int i(0); i<bspl.degree+1; ++i)
    {
        y += bcoef(jmin+i) * vals(i);
    }
    return y;
}

double Spline_1D::eval(double x) const
{
    double values[bspl.degree+1];
    mdspan_1d vals(values, bspl.degree+1);
    
    if (bspl.uniform) return eval_intern<BSplines_uniform>(x, static_cast<const BSplines_uniform&>(bspl), vals);
    else return eval_intern<BSplines_non_uniform>(x, static_cast<const BSplines_non_uniform&>(bspl), vals);
}
    
double Spline_1D::eval_deriv(double x) const
{
    double values[bspl.degree+1];
    mdspan_1d vals(values, bspl.degree+1);
    
    if (bspl.uniform) return eval_deriv_intern<BSplines_uniform>(x, static_cast<const BSplines_uniform&>(bspl), vals);
    else return eval_deriv_intern<BSplines_non_uniform>(x, static_cast<const BSplines_non_uniform&>(bspl), vals);
}

template <class T, typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
void Spline_1D::eval_array_loop(mdspan_1d const& x, mdspan_1d& y) const
{
    const T& l_bspl = static_cast<const T&>(bspl);

    assert( x.extent(0) == y.extent(0) );
    double values[l_bspl.degree+1];
    mdspan_1d vals(values, l_bspl.degree+1);

    for (int i(0); i<x.extent(0); ++i)
    {
        y(i) = eval_intern<T>(x(i), l_bspl, vals);
    }
}

template <class T, typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
void Spline_1D::eval_array_deriv_loop(mdspan_1d const& x, mdspan_1d& y) const
{
    const T& l_bspl = static_cast<const T&>(bspl);

    assert( x.extent(0) == y.extent(0) );
    double values[l_bspl.degree+1];
    mdspan_1d vals(values, l_bspl.degree+1);

    for (int i(0); i<x.extent(0); ++i)
    {
        y(i) = eval_deriv_intern<T>(x(i), l_bspl, vals);
    }
}

void Spline_1D::eval_array(mdspan_1d const& x, mdspan_1d& y) const
{
    if (bspl.uniform) eval_array_loop<BSplines_uniform>(x, y);
    else eval_array_loop<BSplines_non_uniform>(x, y);
}

void Spline_1D::eval_array_deriv(mdspan_1d const& x, mdspan_1d& y) const
{
    if (bspl.uniform) eval_array_deriv_loop<BSplines_uniform>(x, y);
    else eval_array_deriv_loop<BSplines_non_uniform>(x, y);
}

double Spline_1D::integrate() const
{
    double values[bcoef.extent(0)];
    mdspan_1d vals(values, bcoef.extent(0));

    bspl.integrals(vals);

    double y = 0.0;
    for (int i(0); i<bcoef.extent(0); ++i)
    {
        y += bcoef(i) * vals(i);
    }
    return y;
}

