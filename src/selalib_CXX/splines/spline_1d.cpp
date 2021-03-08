#include <cassert>
#include "spline_1d.h"


Spline_1D::Spline_1D(const BSplines* bspl)
    : bcoef_ptr(new double[bspl->degree + bspl->ncells]),
      bcoef(bcoef_ptr, bspl->degree + bspl->ncells),
      bspl(bspl)
{}

Spline_1D::~Spline_1D()
{
    delete[] bcoef_ptr;
}

bool Spline_1D::belongs_to_space(const BSplines* bspline) const
{
    return bspl == bspline;
}

double Spline_1D::eval(double x) const
{
    double values[bspl->degree+1];
    mdspan_1d vals(values, bspl->degree+1);
    int jmin;

    bspl->eval_basis(x, vals, jmin);

    double y = 0.0;
    for (int i(0); i<bspl->degree+1; ++i)
    {
        y += bcoef(jmin+i) * vals(i);
    }
    return y;
}
    
double Spline_1D::eval_deriv(double x) const
{
    double values[bspl->degree+1];
    mdspan_1d derivs(values, bspl->degree+1);
    int jmin;

    bspl->eval_deriv(x, derivs, jmin);

    double y = 0.0;
    for (int i(0); i<bspl->degree+1; ++i)
    {
        y += bcoef(jmin+i) * derivs(i);
    }
    return y;
}

void Spline_1D::eval_array(mdspan_1d const& x, mdspan_1d& y) const
{
    assert( x.extent(0) == y.extent(0) );

    for (int i(0); i<x.extent(0); ++i)
    {
        y(i) = eval(x(i));
    }
}

void Spline_1D::eval_array_deriv(mdspan_1d const& x, mdspan_1d& y) const
{
    assert( x.extent(0) == y.extent(0) );

    for (int i(0); i<x.extent(0); ++i)
    {
        y(i) = eval_deriv(x(i));
    }
}

double Spline_1D::integrate() const
{
    double values[bcoef.extent(0)];
    mdspan_1d vals(values, bcoef.extent(0));

    bspl->integrate(vals);

    double y = 0.0;
    for (int i(0); i<bcoef.extent(0); ++i)
    {
        y += bcoef(i) * vals(i);
    }
    return y;
}

/************************************************************************
 *                    Fortran access functions                          *
 ************************************************************************/
Spline_1D* spline_1d_new(BSplines* bspl)
{
    return new Spline_1D(bspl);
}
void spline_1d_free(Spline_1D* spl)
{
    delete spl;
}
double spline_1d_eval(Spline_1D* spl, double x)
{
    return spl->eval(x);
}
double spline_1d_eval_deriv(Spline_1D* spl, double x)
{
    return spl->eval_deriv(x);
}
void spline_1d_eval_array(Spline_1D* spl, double* x_ptr, int nx, double* y_ptr, int ny)
{
    mdspan_1d x(x_ptr, nx);
    mdspan_1d y(y_ptr, ny);
    spl->eval_array(x,y);
}
void spline_1d_eval_array_deriv(Spline_1D* spl, double* x_ptr, int nx, double* y_ptr, int ny)
{
    mdspan_1d x(x_ptr, nx);
    mdspan_1d y(y_ptr, ny);
    spl->eval_array_deriv(x,y);
}

double spline_1d_integrate(Spline_1D* spl)
{
    return spl->integrate();
}

int spline_1d_get_ncoeffs(Spline_1D* spl)
{
    return spl->bcoef.extent(0);
}

double* spline_1d_get_coeffs(Spline_1D* spl)
{
    return spl->bcoef_ptr;
}

