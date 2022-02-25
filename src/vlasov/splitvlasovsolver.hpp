// SPDX-License-Identifier: MIT

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

    ~SplitVlasovSolver() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_field, double dt) const override;
};
