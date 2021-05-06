#pragma once

#include "block.h"

class IEfieldSolver
{
public:
    virtual void operator()(DBlockViewX& ex, const DBlockViewXVx& fdistribu) const = 0;
};
