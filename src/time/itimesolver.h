#pragma once

#include "blockview.h"

class ITimeSolver
{
public:
    virtual DBlockViewXVx& operator()(DBlockViewXVx& data, double mass_ratio, int steps) const = 0;
};
