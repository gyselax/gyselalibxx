#pragma once

#include "iadvectionx.h"
#include "blockview.h"
#include "ivlasovsolver.h"

class VlasovSolver : public IVlasovSolver
{
    const IAdvectionX& m_advec_x;

    const IAdvectionX& m_advec_vx;

public:
    VlasovSolver(const IAdvectionX& advec_x, const IAdvectionX& m_advec_vx);

    DBlockViewXVx& operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt) const override;
};
