#pragma once

#include <geometry.hpp>

class ITimeSolver
{
public:
    virtual void operator()(DSpanSpXVx allfdistribu, int steps) const = 0;
};
