#include <cassert>
#include "bsplines.h"
#include "bsplines_non_uniform.h"
#include "bsplines_uniform.h"

BSplines::BSplines(int degree, bool periodic, bool uniform, int ncells, int nbasis, double xmin, double xmax, bool radial)
    : degree(degree), periodic(periodic), uniform(uniform), radial(radial), ncells(ncells), nbasis(nbasis), xmin(xmin), xmax(xmax)
{}

BSplines* BSplines::new_bsplines(int degree, bool periodic, std::vector<double> const& breaks)
{
    return new BSplines_non_uniform(degree, periodic, breaks);
}
    
BSplines* BSplines::new_bsplines(int degree, bool periodic, double xmin, double xmax, int ncells)
{
    return new BSplines_uniform(degree, periodic, xmin, xmax, ncells);
}

BSplines* get_new_bspline_uniform(int degree, bool periodic, double xmin, double xmax, int ncells)
{
    return BSplines::new_bsplines(degree, periodic, xmin, xmax, ncells);
}

BSplines* get_new_bspline_non_uniform(int degree, bool periodic, double xmin, double xmax, int ncells, double* breaks_ptr, int nbreaks)
{
    std::vector<double> breaks;
    assert( nbreaks == ncells+1 );
    breaks.reserve(ncells+1);
    for (int i(0); i<(ncells+1); ++i)
    {
        breaks.push_back(breaks_ptr[i]);
    }
    return BSplines::new_bsplines(degree, periodic, breaks);
}

void bsplines_free(BSplines* bspl)
{
    delete bspl;
}

void bsplines_eval_basis(BSplines* bspl, double x, int nvals, double* vals, int* jmin)
{
    mdspan_1d values(vals, nvals);
    bspl->eval_basis(x, values, *jmin);
}

void bsplines_eval_deriv(BSplines* bspl, double x, int nvals, double* vals, int* jmin)
{
    mdspan_1d values(vals, nvals);
    bspl->eval_deriv(x, values, *jmin);
}

void bsplines_eval_basis_and_n_derivs(BSplines* bspl, double x, int n, int nvals1, int nvals2, double* vals, int* jmin)
{
    mdspan_2d values(vals, nvals2, nvals1);
    bspl->eval_basis_and_n_derivs(x, n, values, *jmin);
}

void bsplines_integrals(BSplines* bspl, int nvals, double* int_vals)
{
    mdspan_1d values(int_vals, nvals);
    bspl->integrate(values);
}

double bsplines_get_knot(BSplines* bspl, int idx)
{
    return bspl->get_knot(idx);
}
