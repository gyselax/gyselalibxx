#include "sll/deprecated/bsplines.hpp"

namespace deprecated {

BSplines::BSplines(
        int degree,
        bool periodic,
        int ncells,
        int nbasis,
        double xmin,
        double xmax,
        bool radial)
    : m_degree(degree)
    , m_radial(radial)
    , m_nbasis(nbasis)
    , m_ncells(ncells)
    , m_periodic(periodic)
    , m_xmin(xmin)
    , m_xmax(xmax)
    , m_length(xmax - xmin)
{
}

} // namespace deprecated
