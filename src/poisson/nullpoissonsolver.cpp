#include "nullpoissonsolver.hpp"

DSpanX NullPoissonSolver::operator()(DSpanX const electostatic_potential, DViewSpXVx) const
{
    return electostatic_potential;
}
