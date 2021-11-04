#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class IAdvectionVx
{
public:
    virtual ~IAdvectionVx() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const = 0;
};
