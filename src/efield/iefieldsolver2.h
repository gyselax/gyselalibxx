#pragma once

#include "blockview.h"

class IEfieldSolver2
{
public:
    virtual DBlockViewX& operator()(DBlockViewX& ex, const DBlockViewXVx& fdistribu) const = 0;
};
