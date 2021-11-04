#pragma once

#include <geometry.hpp>

#include "ipoissonsolver.hpp"

class NullPoissonSolver : public IPoissonSolver
{
public:
    NullPoissonSolver() = default;

    ~NullPoissonSolver() override = default;

    DSpanX operator()(DSpanX electic_potential, DViewSpXVx allfdistribu) const override;
};
