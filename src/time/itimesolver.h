#pragma once

#include "blockview.h"

class ITimeSolver
{
public:
    virtual DBlockSpanXVx operator()(DBlockSpanXVx data, double mass_ratio, int steps) const = 0;
};
