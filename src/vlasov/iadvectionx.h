#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IAdvectionX
{
public:
    virtual DSpanXVx operator()(DSpanXVx fdistribu, double sqrt_me_on_mspecies, double dt)
            const = 0;
};
