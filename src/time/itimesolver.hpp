// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, int steps) const = 0;
};
