#pragma once

#include <ddc/BlockSpan>

#include <geometry.hpp>

class IPoissonSolver
{
public:
    virtual DSpanX operator()(DSpanX electric_potential, DViewSpXVx allfdistribu) const = 0;
};
