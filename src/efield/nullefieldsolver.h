#pragma once

#include "iefieldsolver.h"

class NullEfieldSolver : public IEfieldSolver
{
public:
    DBlockViewX& operator()(DBlockViewX& ex, const DBlockViewXVx& fdistribu) const override;
};
