#pragma once

#include "fdistribu.h"

class ITimeSolver
{
public:
    virtual void operator()(DistributionFunction& fdistribu, double mass_ratio, int steps)
            const = 0;
};
