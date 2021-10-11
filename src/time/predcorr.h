#pragma once

#include <geometry.h>

#include "fdistribu.h"
#include "itimesolver.h"

class IEfieldSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
    IVlasovSolver const& m_vlasov_solver;

    IEfieldSolver const& m_efield_solver;

    double const m_dt;

    double const m_time_diag;

public:
    PredCorr(
            const IVlasovSolver& vlasov_solver,
            const IEfieldSolver& efield_solver,
            double dt,
            double time_diag);

    void operator()(DistributionFunction& fdistribu, double electron_mass, int steps)
            const override;
};
