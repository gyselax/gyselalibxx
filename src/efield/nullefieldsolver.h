#pragma once

#include "iefieldsolver.h"

class NullEfieldSolver : public IEfieldSolver
{
public:
    DBlockSpanX operator()(DBlockSpanX ex, DBlockViewXVx fdistribu) const override;
};
