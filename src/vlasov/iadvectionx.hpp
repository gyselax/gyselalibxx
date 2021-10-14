#pragma once

#include <ddc/BlockSpan>

#include <geometry.hpp>

class IAdvectionX
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const = 0;
};
