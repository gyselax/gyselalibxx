#pragma once

#include <geometry.h>

#include "itimesolver.h"

class IEfieldSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
    IVlasovSolver const& m_vlasov;

    IEfieldSolver const& m_efield;

    RLengthT const m_dt;

public:
    PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver& efield, RLengthT dt);

    DSpanXVx operator()(DSpanXVx data, double mass_ratio, int steps) const override;
};
