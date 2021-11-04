#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class IPoissonSolver
{
public:
    virtual ~IPoissonSolver() = default;

    virtual DSpanX operator()(DSpanX electric_potential, DViewSpXVx allfdistribu) const = 0;
};
