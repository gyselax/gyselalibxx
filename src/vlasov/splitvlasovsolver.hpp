#pragma once

#include <geometry.hpp>

#include "ivlasovsolver.hpp"

class IAdvectionVx;
class IAdvectionX;

class SplitVlasovSolver : public IVlasovSolver
{
    IAdvectionX const& m_advec_x;

    IAdvectionVx const& m_advec_vx;

public:
    SplitVlasovSolver(IAdvectionX const& advec_x, IAdvectionVx const& m_advec_vx);

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX efield, double dt) const override;
};
