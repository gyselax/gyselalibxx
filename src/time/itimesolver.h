#pragma once

#include <geometry.h>

class ITimeSolver
{
public:
    virtual void operator()(DSpanSpXVx allfdistribu, int steps) const = 0;
};
