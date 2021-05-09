#pragma once

#include "iefieldsolver2.h"
#include "itimesolver2.h"
#include "ivlasovsolver2.h"

class PredCorr2 : public ITimeSolver2
{
    IVlasovSolver2 const& m_vlasov;

    IEfieldSolver2 const& m_efield;

    RLengthT const m_dt;

public:
    PredCorr2(const IVlasovSolver2& vlasov, const IEfieldSolver2& efield, RLengthT dt);

    DBlockViewXVx& operator()(DBlockViewXVx& data, double mass_ratio, int steps) const override;
};
