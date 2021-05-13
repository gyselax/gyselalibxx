#pragma once

#include "blockview.h"
#include "iadvectionx.h"
#include "ivlasovsolver.h"

class SplitVlasovSolver : public IVlasovSolver
{
    const IAdvectionX& m_advec_x;

    const IAdvectionX& m_advec_vx;

public:
    SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionX& m_advec_vx);

    DBlockViewXVx& operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt)
            const override;
};
