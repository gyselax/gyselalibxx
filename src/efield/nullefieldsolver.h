#pragma once

#include "iefieldsolver.h"

class NullEfieldSolver : public IEfieldSolver
{
public:
    DBlockViewX operator()(DBlockViewX ex, DBlockCViewXVx fdistribu) const override;
};
