#pragma once

#include "advection1d.h"
#include "block.h"
#include "ivlasovsolver.h"

class VlasovSolver : public IVlasovSolver
{
    const Advection1D& m_advec_x;

    const Advection1D& m_advec_vx;

public:
    VlasovSolver(const Advection1D& advec_x, const Advection1D& m_advec_vx);

    void operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt) const override;
};
