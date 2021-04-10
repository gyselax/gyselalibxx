#pragma once

#include "block.h"

class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    virtual void operator()(DBlock2D& data, double mass_ratio) const = 0;
};
