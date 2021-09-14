#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IAdvectionVx
{
public:
    virtual DSpanXVx operator()(DSpanXVx fdistribu, double mass_ratio, double dt) const = 0;
};
