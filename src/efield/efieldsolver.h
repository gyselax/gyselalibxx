#pragma once

#include "iefieldsolver.h"

class EfieldSolver : public IEfieldSolver
{
public:
    void operator()(DBlock1D ex, const DBlock2D fdistribu) const override;
};
