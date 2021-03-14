#pragma once

#include "block.h"

class EfieldSolver
{
public:
    virtual ~EfieldSolver() = default;

    virtual void operator()(DBlock1D ex, const DBlock2D fdistribu) const = 0;
};
