#pragma once

#include <geometry.hpp>

#include "ivlasovsolver.hpp"

class IAdvectionVx;
class IAdvectionX;

class SplitVlasovSolver : public IVlasovSolver
{
    const IAdvectionX& m_advec_x;

    const IAdvectionVx& m_advec_vx;

public:
    SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionVx& m_advec_vx);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX efield, double dt) const override;
};
