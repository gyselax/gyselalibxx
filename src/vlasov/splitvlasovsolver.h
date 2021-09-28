#pragma once

#include <geometry.h>

#include "ivlasovsolver.h"

class IAdvectionVx;
class IAdvectionX;

class SplitVlasovSolver : public IVlasovSolver
{
    const IAdvectionX& m_advec_x;

    const IAdvectionVx& m_advec_vx;

public:
    SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionVx& m_advec_vx);

    DSpanXVx operator()(DSpanXVx fdistribu, DViewX efield, double sqrt_me_on_mspecies, double dt)
            const override;
};
