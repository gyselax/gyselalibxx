// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <geometry.hpp>

class IPoissonSolver
{
public:
    virtual ~IPoissonSolver() = default;

    virtual void operator()(
            DSpanRP electrostatic_potential,
            DSpanRP electric_field,
            DViewRP allfdistribu) const = 0;
};
