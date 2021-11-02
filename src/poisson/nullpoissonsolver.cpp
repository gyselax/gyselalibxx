#include "nullpoissonsolver.hpp"

DSpanX NullPoissonSolver::operator()(
        DSpanX const electostatic_potential,
        DViewSpXVx const allfdistribu) const
{
    return electostatic_potential;
}
