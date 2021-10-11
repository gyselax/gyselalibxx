#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IEfieldSolver
{
public:
    virtual DSpanX operator()(DSpanX efield, DViewXVx fdistribu) const = 0;
};
