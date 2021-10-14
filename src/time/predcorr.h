#pragma once

#include <geometry.h>

#include "itimesolver.h"

class IPoissonSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
    IVlasovSolver const& m_vlasov_solver;

    IPoissonSolver const& m_poisson_solver;

    double const m_dt;

public:
    PredCorr(const IVlasovSolver& vlasov_solver, const IPoissonSolver& poisson_solver, double dt);

    void operator()(DSpanSpXVx allfdistribu, int steps) const override;
};
