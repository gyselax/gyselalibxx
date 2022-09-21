// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "itimesolver.hpp"

class IPoissonSolver;
class IBoltzmannSolver;

class PredCorr : public ITimeSolver
{
    IBoltzmannSolver const& m_boltzmann_solver;

    IPoissonSolver const& m_poisson_solver;

    double const m_dt;

public:
    PredCorr(
            IBoltzmannSolver const& boltzmann_solver,
            IPoissonSolver const& poisson_solver,
            double dt);

    ~PredCorr() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, int steps) const override;
};
