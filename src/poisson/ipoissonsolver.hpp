#pragma once

#include <ddc/ChunkSpan>

#include <geometry.hpp>

class SpeciesInformation;

class IPoissonSolver
{
public:
    virtual ~IPoissonSolver() = default;

    virtual void operator()(
            DSpanX electrostatic_potential,
            DSpanX electric_field,
            DViewSpXVx allfdistribu) const = 0;
};
