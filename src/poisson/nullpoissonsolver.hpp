#pragma once

#include <geometry.hpp>

#include "ipoissonsolver.hpp"

class NullPoissonSolver : public IPoissonSolver
{
public:
    DSpanX operator()(DSpanX electic_potential, DViewSpXVx allfdistribu) const override;
};
