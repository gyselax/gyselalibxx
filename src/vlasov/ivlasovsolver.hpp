#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class IVlasovSolver
{
public:
    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX efield, double dt) const = 0;
};
