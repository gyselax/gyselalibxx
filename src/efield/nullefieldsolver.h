#pragma once

#include "iefieldsolver2.h"

class NullEfieldSolver : public IEfieldSolver2
{
public:
    DBlockViewX& operator()(DBlockViewX& ex, const DBlockViewXVx& fdistribu) const override;
};
