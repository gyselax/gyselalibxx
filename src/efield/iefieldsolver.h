#pragma once

#include "blockview.h"

class IEfieldSolver
{
public:
    virtual DBlockSpanX operator()(DBlockSpanX ex, DBlockViewXVx fdistribu) const = 0;
};
