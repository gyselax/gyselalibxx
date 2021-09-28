#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IAdvectionVx
{
public:
    virtual DSpanXVx operator()(
            DSpanXVx fdistribu,
            DViewX efield,
            double sqrt_me_on_mspecies,
            double dt) const = 0;
};
