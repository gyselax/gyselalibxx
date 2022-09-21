// SPDX-License-Identifier: MIT

#pragma once

#include <utility>
#include <vector>

#include <geometry.hpp>
#include <isource.hpp>

#include "iboltzmannsolver.hpp"

class SplitSourceEnvironmentSolver : public IBoltzmannSolver
{
    IBoltzmannSolver const& m_boltzmann_solver;

    std::vector<std::reference_wrapper<ISource const>> m_rhs;

public:
    SplitSourceEnvironmentSolver(
            IBoltzmannSolver const& vlasov_solver,
            std::vector<std::reference_wrapper<ISource const>> rhs);

    ~SplitSourceEnvironmentSolver() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_field, double dt) const override;
};
