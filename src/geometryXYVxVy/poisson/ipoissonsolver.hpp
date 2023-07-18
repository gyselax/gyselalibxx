// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

class IPoissonSolver
{
public:
    virtual ~IPoissonSolver() = default;

    virtual void operator()(
            DSpanXY electrostatic_potential,
            DSpanXY electric_field_x,
            DSpanXY electric_field_y,
            DViewSpXYVxVy allfdistribu) const = 0;
};
