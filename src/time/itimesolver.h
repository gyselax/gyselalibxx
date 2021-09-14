#pragma once

#include "block_span.h"

class ITimeSolver
{
public:
    virtual DSpanXVx operator()(DSpanXVx data, double mass_ratio, int steps) const = 0;
};
