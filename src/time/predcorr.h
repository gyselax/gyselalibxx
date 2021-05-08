#pragma once

#include "iefieldsolver2.h"
#include "itimesolver.h"
#include "ivlasovsolver2.h"

class PredCorr : public ITimeSolver
{
    const IVlasovSolver2& m_vlasov;

    const IEfieldSolver2& m_efield;

    const MDomain<Dim::T>& m_time;

public:
    PredCorr(const IVlasovSolver2& vlasov, const IEfieldSolver2& efield, const MDomain<Dim::T>& time);

    DBlockViewXVx& operator()(DBlockViewXVx& data, double mass_ratio) const override;
};
