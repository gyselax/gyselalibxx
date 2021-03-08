#include "spline_2d.h"
#include <cassert>
#include <iostream>


Spline_2D::Spline_2D(BSplines* bspl1, BSplines* bspl2)
    : bcoef_ptr(new double[(bspl1->degree + bspl1->ncells)*(bspl2->degree + bspl2->ncells)]),
      bcoef(bcoef_ptr, bspl1->degree + bspl1->ncells, bspl2->degree + bspl2->ncells),
      bspl1(bspl1),
      bspl2(bspl2)
{}

Spline_2D::~Spline_2D()
{
    delete[] bcoef_ptr;
}

bool Spline_2D::belongs_to_space(const BSplines* bspline1, const BSplines* bspline2) const
{
    return (bspl1 == bspline1 && bspl2 == bspline2);
}

double Spline_2D::eval(const double x1, const double x2) const
{
    double values1[bspl1->degree+1];
    double values2[bspl2->degree+1];
    mdspan_1d vals1(values1, bspl1->degree+1);
    mdspan_1d vals2(values2, bspl2->degree+1);
    int jmin1, jmin2;

    bspl1->eval_basis(x1, vals1, jmin1);
    bspl2->eval_basis(x2, vals2, jmin2);

    double y = 0.0;
    for (int i(0); i<bspl1->degree+1; ++i)
    {
        for (int j(0); j<bspl2->degree+1; ++j)
        {
            y += bcoef(jmin1 + i, jmin2 + j) * vals1(i) * vals2(j);
        }
    }
    return y;
}

double Spline_2D::eval_deriv(const double x1, const double x2, const bool deriv1, const bool deriv2) const
{
    double values1[bspl1->degree+1];
    double values2[bspl2->degree+1];
    mdspan_1d vals1(values1, bspl1->degree+1);
    mdspan_1d vals2(values2, bspl2->degree+1);

    int jmin1, jmin2;

    if (deriv1)
    {
        bspl1->eval_deriv(x1, vals1, jmin1);
    }
    else
    {
        bspl1->eval_basis(x1, vals1, jmin1);
    }
    if (deriv2)
    {
        bspl2->eval_deriv(x2, vals2, jmin2);
    }
    else
    {
        bspl2->eval_basis(x2, vals2, jmin2);
    }

    double y = 0.0;
    for (int i(0); i<bspl1->degree+1; ++i)
    {
        for (int j(0); j<bspl2->degree+1; ++j)
        {
            y += bcoef(jmin1 + i, jmin2 + j) * vals1(i) * vals2(j);
        }
    }
    return y;
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

void Spline_2D::eval_array_deriv(mdspan_2d const& x1, mdspan_2d const& x2, mdspan_2d& y, const bool deriv1, const bool deriv2) const
{
    assert( x1.extent(0) == y.extent(0) );
    assert( x1.extent(1) == y.extent(1) );
    assert( x2.extent(0) == y.extent(0) );
    assert( x2.extent(1) == y.extent(1) );

    for (int i(0); i<x1.extent(0); ++i)
    {
        for (int j(0); j<x1.extent(1); ++j)
        {
            y(i,j) = eval_deriv(x1(i,j), x2(i,j), deriv1, deriv2);
        }
    }
}

void Spline_2D::integrate_dim(mdspan_1d& y, const int dim) const
{
    assert( dim >=0 and dim < 2 );
    assert( y.extent(0) == bcoef.extent(1 - dim) );

    const BSplines* bspline ((dim == 0) ? bspl1 : bspl2);

    double values[bcoef.extent(dim)];
    mdspan_1d vals(values, bcoef.extent(dim));

    bspline->integrate(vals);

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

    bspl1->integrate(vals);

    double y = 0.0;
    for (int i(0); i<bcoef.extent(0); ++i)
    {
        y += int_vals(i) * vals(i);
    }
    return y;
}

/************************************************************************
 *                    Fortran access functions                          *
 ************************************************************************/
Spline_2D* spline_2d_new(BSplines* bspl1, BSplines* bspl2)
{
    return new Spline_2D(bspl2, bspl1);
}
void spline_2d_free(Spline_2D* spl)
{
    delete spl;
}
double spline_2d_eval(Spline_2D* spl, double* x1, double* x2)
{
    return spl->eval(*x2, *x1);
}
double spline_2d_eval_deriv(Spline_2D* spl, double* x1, double* x2, bool deriv_x1, bool deriv_x2)
{
    return spl->eval_deriv(*x2, *x1, deriv_x2, deriv_x1);
}
void spline_2d_eval_array(Spline_2D* spl,
        double* x1_ptr, int nx1_1, int nx1_2,
        double* x2_ptr, int nx2_1, int nx2_2,
        double* y_ptr,  int ny1,   int ny2)
{
    mdspan_2d x1(x1_ptr, nx1_2, nx1_1);
    mdspan_2d x2(x2_ptr, nx2_2, nx2_1);
    mdspan_2d y(y_ptr, ny2, ny1);
    spl->eval_array(x2,x1,y);
}
void spline_2d_eval_array_deriv(Spline_2D* spl,
        double* x1_ptr, int nx1_1, int nx1_2,
        double* x2_ptr, int nx2_1, int nx2_2,
        bool deriv_x1, bool deriv_x2,
        double* y_ptr, int ny1, int ny2)
{
    mdspan_2d x1(x1_ptr, nx1_2, nx1_1);
    mdspan_2d x2(x2_ptr, nx2_2, nx2_1);
    mdspan_2d y(y_ptr, ny2, ny1);
    spl->eval_array_deriv(x2,x1,y,deriv_x2, deriv_x1);
}

void spline_2d_integrate_dim(Spline_2D* spl, double* y_ptr, int ny, const int dim)
{
    mdspan_1d y(y_ptr, ny);
    spl->integrate_dim(y, dim);
}
double spline_2d_integrate(Spline_2D* spl)
{
    return spl->integrate();
}
