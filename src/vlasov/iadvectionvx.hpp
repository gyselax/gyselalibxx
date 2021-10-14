#pragma once

#include <ddc/BlockSpan>

#include <geometry.hpp>

class IAdvectionVx
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_potential, double dt)
            const = 0;
};
