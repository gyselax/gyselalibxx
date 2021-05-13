#pragma once

#include "blockview.h"

class IEfieldSolver
{
public:
    virtual DBlockViewX operator()(DBlockViewX ex, DBlockCViewXVx fdistribu) const = 0;
};
