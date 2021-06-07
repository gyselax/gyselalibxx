#pragma once

#include "blockview.h"
#include "geometry.h"

class IEfieldSolver
{
public:
    virtual DBlockSpanX operator()(DBlockSpanX ex, DBlockViewXVx fdistribu) const = 0;
};
