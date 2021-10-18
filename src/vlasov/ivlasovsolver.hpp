#pragma once

#include <ddc/ChunkSpan>

class IVlasovSolver
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX efield, double dt) const = 0;
};
