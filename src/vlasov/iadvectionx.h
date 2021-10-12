#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IAdvectionX
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx fdistribu, double dt) const = 0;
};
