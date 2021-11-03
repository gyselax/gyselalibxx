#pragma once

#include <geometry.hpp>

class ITimeSolver
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, int steps) const = 0;
};
