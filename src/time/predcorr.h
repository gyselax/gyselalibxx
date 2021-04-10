#pragma once

#include "iefieldsolver.h"
#include "itimesolver.h"
#include "ivlasovsolver.h"

class PredCorr : public ITimeSolver
{
    const IVlasovSolver& m_vlasov;

    const IEfieldSolver& m_efield;

    const MDomain1D& m_time;

public:
    PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver& efield, const MDomain1D& time);

    void operator()(DBlock2D& data, double mass_ratio) const override;
};
