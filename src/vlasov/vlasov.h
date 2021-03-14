#pragma once

#include "advection1d.h"
#include "block.h"

class Vlasov
{
    const Advection1D& m_advec_x;

    const Advection1D& m_advec_vx;

public:
    Vlasov(const Advection1D& advec_x, const Advection1D& m_advec_vx);

    void operator()(DBlock2D& cur, double dt) const;
};
