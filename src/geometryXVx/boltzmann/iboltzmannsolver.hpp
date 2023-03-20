// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IBoltzmannSolver
{
public:
    virtual ~IBoltzmannSolver() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX efield, double dt) const = 0;
};
