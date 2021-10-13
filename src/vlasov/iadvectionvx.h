#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IAdvectionVx
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx fdistribu, DViewX electric_potential, double dt)
            const = 0;
};
