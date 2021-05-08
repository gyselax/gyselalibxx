#pragma once

#include "iefieldsolver2.h"
#include "itimesolver2.h"
#include "ivlasovsolver2.h"

class PredCorr2 : public ITimeSolver2
{
    const IVlasovSolver2& m_vlasov;

    const IEfieldSolver2& m_efield;

    const MDomain<Dim::T>& m_time;

public:
    PredCorr2(
            const IVlasovSolver2& vlasov,
            const IEfieldSolver2& efield,
            const MDomain<Dim::T>& time);

    DBlockViewXVx& operator()(DBlockViewXVx& data, double mass_ratio) const override;
};
