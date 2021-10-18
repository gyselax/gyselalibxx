#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class IAdvectionX
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const = 0;
};
