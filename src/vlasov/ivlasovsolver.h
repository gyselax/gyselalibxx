#pragma once

#include <ddc/BlockSpan>

class IVlasovSolver
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx fdistribu, DViewX efield, double dt) const = 0;
};
