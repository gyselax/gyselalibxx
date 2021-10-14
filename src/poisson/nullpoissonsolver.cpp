#include "nullpoissonsolver.hpp"

DSpanX NullPoissonSolver::operator()(DSpanX electostatic_potential, DViewSpXVx allfdistribu) const
{
    return electostatic_potential;
}
