#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class IVlasovSolver
{
public:
    virtual ~IVlasovSolver() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX efield, double dt) const = 0;
};
