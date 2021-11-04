#pragma once

#include <geometry.hpp>

#include "itimesolver.hpp"

class IPoissonSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
    IVlasovSolver const& m_vlasov_solver;

    IPoissonSolver const& m_poisson_solver;

    double const m_dt;

public:
    PredCorr(IVlasovSolver const& vlasov_solver, IPoissonSolver const& poisson_solver, double dt);

    ~PredCorr() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, int steps) const override;
};
