#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IAdvectionVx
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx fdistribu, DViewX efield, double dt) const = 0;
};
