#pragma once

#include "blockview.h"
#include "geometry.h"

class IAdvectionX
{
public:
    virtual DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt)
            const = 0;
};
