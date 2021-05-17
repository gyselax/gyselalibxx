#pragma once

#include "iefieldsolver.h"
#include "itimesolver.h"
#include "ivlasovsolver.h"

class PredCorr : public ITimeSolver
{
    IVlasovSolver const& m_vlasov;

    IEfieldSolver const& m_efield;

    RLengthT const m_dt;

public:
    PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver& efield, RLengthT dt);

    DBlockSpanXVx operator()(DBlockSpanXVx data, double mass_ratio, int steps) const override;
};
