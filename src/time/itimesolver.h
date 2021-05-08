#pragma once

#include "blockview.h"

class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    virtual DBlockViewXVx& operator()(DBlockViewXVx& data, double mass_ratio) const = 0;
};
