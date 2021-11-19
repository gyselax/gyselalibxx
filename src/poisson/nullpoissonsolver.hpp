#pragma once

#include <geometry.hpp>

#include "ipoissonsolver.hpp"

class NullPoissonSolver : public IPoissonSolver
{
public:
    NullPoissonSolver() = default;

    ~NullPoissonSolver() override = default;

    void operator()(DSpanX electrostatic_potential, DSpanX electric_field, DViewSpXVx allfdistribu)
            const override;
};
