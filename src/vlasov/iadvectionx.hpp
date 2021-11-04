#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class IAdvectionX
{
public:
    virtual ~IAdvectionX() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const = 0;
};
