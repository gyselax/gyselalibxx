// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt, int steps = 1) const = 0;
};
