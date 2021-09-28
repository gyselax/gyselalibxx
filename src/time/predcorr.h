#pragma once

#include <geometry.h>

#include "fdistribu.h"
#include "itimesolver.h"

class IEfieldSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
    IVlasovSolver const& m_vlasov;

    IEfieldSolver const& m_efield;

    double const m_dt;

public:
    PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver& efield, double dt);

    void operator()(DistributionFunction& fdistribu, double mass_ratio, int steps) const override;
};
