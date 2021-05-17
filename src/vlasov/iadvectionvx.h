#pragma once

#include "blockview.h"

class IAdvectionVx
{
public:
    virtual DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt)
            const = 0;
};
