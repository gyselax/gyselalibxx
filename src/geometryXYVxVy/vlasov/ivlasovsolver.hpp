// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IVlasovSolver
{
public:
    virtual ~IVlasovSolver() = default;

    virtual DSpanSpXYVxVy operator()(
            DSpanSpXYVxVy allfdistribu,
            DViewXY efield_x,
            DViewXY efield_y,
            double dt) const = 0;
};
