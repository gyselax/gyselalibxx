#pragma once

#include <geometry.h>

#include "itimesolver.h"

class IEfieldSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
    IVlasovSolver const& m_vlasov_solver;

    IEfieldSolver const& m_efield_solver;

    double const m_dt;

public:
    PredCorr(const IVlasovSolver& vlasov_solver, const IEfieldSolver& efield_solver, double dt);

    void operator()(DSpanSpXVx fdistribu, int steps) const override;
};
