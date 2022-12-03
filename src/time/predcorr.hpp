// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "itimesolver.hpp"

class IPoissonSolver;
class IBoltzmannSolver;

class PredCorr : public ITimeSolver
{
private:
    IBoltzmannSolver const& m_boltzmann_solver;

    IPoissonSolver const& m_poisson_solver;

public:
    PredCorr(IBoltzmannSolver const& boltzmann_solver, IPoissonSolver const& poisson_solver);

    ~PredCorr() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt, int steps = 1) const override;
};
