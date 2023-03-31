// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ipoissonsolver.hpp"

class NullPoissonSolver : public IPoissonSolver
{
public:
    NullPoissonSolver() = default;

    ~NullPoissonSolver() override = default;

    void operator()(DSpanRP electrostatic_potential, DSpanRP electric_field, DViewRP allfdistribu)
            const override;
};
