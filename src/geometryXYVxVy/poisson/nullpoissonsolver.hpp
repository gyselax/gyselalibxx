// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ipoissonsolver.hpp"

class NullPoissonSolver : public IPoissonSolver
{
public:
    NullPoissonSolver() = default;

    ~NullPoissonSolver() override = default;

    void operator()(
            DSpanXY electrostatic_potential,
            DSpanXY electric_field_x,
            DSpanXY electric_field_y,
            DViewSpXYVxVy allfdistribu) const override;
};
