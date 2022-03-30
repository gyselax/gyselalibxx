// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

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
