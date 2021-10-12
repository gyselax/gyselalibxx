#pragma once

#include <geometry.h>

class ITimeSolver
{
public:
    virtual void operator()(DSpanSpXVx fdistribu, int steps) const = 0;
};
