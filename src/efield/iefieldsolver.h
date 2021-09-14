#pragma once

#include <ddc/BlockSpan>

#include <geometry.h>

class IEfieldSolver
{
public:
    virtual DSpanX operator()(DSpanX ex, DViewXVx fdistribu) const = 0;
};
