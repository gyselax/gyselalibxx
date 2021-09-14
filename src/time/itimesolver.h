#pragma once

#include <ddc/BlockSpan>

class ITimeSolver
{
public:
    virtual DSpanXVx operator()(DSpanXVx data, double mass_ratio, int steps) const = 0;
};
