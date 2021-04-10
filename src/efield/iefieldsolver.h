#pragma once

#include "block.h"

class IEfieldSolver
{
public:
    virtual ~IEfieldSolver() = default;

    virtual void operator()(DBlock1D ex, const DBlock2D fdistribu) const = 0;
};
