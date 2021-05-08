#pragma once

#include "iefieldsolver2.h"
#include "itimesolver.h"
#include "ivlasovsolver.h"

class PredCorr : public ITimeSolver
{
    const IVlasovSolver& m_vlasov;

    const IEfieldSolver2& m_efield;

    const MDomain<Dim::T>& m_time;

public:
    PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver2& efield, const MDomain<Dim::T>& time);

    DBlockViewXVx& operator()(DBlockViewXVx& data, double mass_ratio) const override;
};
