// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "itimesolver.hpp"

class IPoissonSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
private:
    IVlasovSolver const& m_vlasov_solver;

    IPoissonSolver const& m_poisson_solver;

public:
    PredCorr(IVlasovSolver const& vlasov_solver, IPoissonSolver const& poisson_solver);

    ~PredCorr() override = default;

    DSpanSpXYVxVy operator()(DSpanSpXYVxVy allfdistribu, double dt, int steps = 1) const override;
};
