#pragma once

#include "blockview.h"

class IEfieldSolver
{
public:
    virtual DBlockViewX& operator()(DBlockViewX& ex, const DBlockViewXVx& fdistribu) const = 0;
};
