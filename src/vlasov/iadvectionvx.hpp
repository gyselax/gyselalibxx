// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class IAdvectionVx
{
public:
    virtual ~IAdvectionVx() = default;

    virtual DSpanSpXVx operator()(
            DSpanSpXVx allfdistribu,
            DViewX electrostatic_potential,
            double dt) const = 0;
};
