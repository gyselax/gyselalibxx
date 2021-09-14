#pragma once

#include "geometry.h"
#include "iefieldsolver.h"

class NullEfieldSolver : public IEfieldSolver
{
public:
    DSpanX operator()(DSpanX ex, DViewXVx fdistribu) const override;
};
