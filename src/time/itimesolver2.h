#pragma once

#include "blockview.h"

class ITimeSolver2
{
public:
    virtual ~ITimeSolver2() = default;

    virtual DBlockViewXVx& operator()(DBlockViewXVx& data, double mass_ratio) const = 0;
};
