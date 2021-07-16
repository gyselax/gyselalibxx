#pragma once

#include "geometry.h"
#include "ivlasovsolver.h"

class IAdvectionVx;
class IAdvectionX;

class SplitVlasovSolver : public IVlasovSolver
{
    const IAdvectionX& m_advec_x;

    const IAdvectionVx& m_advec_vx;

public:
    SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionVx& m_advec_vx);

    DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt) const override;
};
