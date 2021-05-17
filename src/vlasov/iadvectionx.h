#pragma once

#include "blockview.h"

class IAdvectionX
{
public:
    virtual DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt)
            const = 0;
};
