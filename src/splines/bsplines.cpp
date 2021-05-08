#include <cassert>

#include "bsplines.h"

BSplines::BSplines(
        int degree,
        bool periodic,
        bool uniform,
        int ncells,
        int nbasis,
        double xmin,
        double xmax,
        bool radial)
    : degree(degree)
    , periodic(periodic)
    , uniform(uniform)
    , radial(radial)
    , ncells(ncells)
    , nbasis(nbasis)
    , xmin(xmin)
    , xmax(xmax)
    , length(xmax - xmin)
{
}

BSplines* get_new_bspline_uniform(int degree, bool periodic, double xmin, double xmax, int ncells)
{
    return BSplines::new_bsplines(degree, periodic, xmin, xmax, ncells);
}

BSplines* get_new_bspline_non_uniform(
        int degree,
        bool periodic,
        double xmin,
        double xmax,
        int ncells,
        double* breaks_ptr,
        int nbreaks)
{
    std::vector<double> breaks;
    assert(nbreaks == ncells + 1);
    breaks.reserve(ncells + 1);
    for (int i(0); i < (ncells + 1); ++i) {
        breaks.push_back(breaks_ptr[i]);
    }
    return BSplines::new_bsplines(degree, periodic, breaks);
}
