#pragma once

#include <ddc/block_span.h>

#include <geometry.h>

class IEfieldSolver
{
public:
    virtual DSpanX operator()(DSpanX ex, DViewXVx fdistribu) const = 0;
};
