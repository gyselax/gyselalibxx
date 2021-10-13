#include "nullpoissonsolver.hpp"

DSpanX NullPoissonSolver::operator()(DSpanX electostatic_potential, DViewSpXVx fdistribu) const
{
    return electostatic_potential;
}
