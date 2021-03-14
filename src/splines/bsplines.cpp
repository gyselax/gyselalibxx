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
