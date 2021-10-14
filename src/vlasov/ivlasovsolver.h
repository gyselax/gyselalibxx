#pragma once

#include <ddc/BlockSpan>

class IVlasovSolver
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX efield, double dt) const = 0;
};
