#pragma once

#include "block.h"
#include "vlasov.h"

class PredCorr
{
    const Vlasov& m_vlasov;

    const EfieldSolver& m_efield;

public:
    PredCorr(const Vlasov& vlasov, const EfieldSolver& efield);

    void operator()(DBlock2D& data, const MDomain1D& time) const;
};
