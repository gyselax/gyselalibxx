#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IAdvectionX
{
public:
    virtual DSpanXVx operator()(DSpanXVx fdistribu, double mass_ratio, double dt) const = 0;
};
