#pragma once

#include "block.h"

class IAdvectionX
{
public:
    virtual DBlockViewXVx& operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt) const = 0;
    
};


