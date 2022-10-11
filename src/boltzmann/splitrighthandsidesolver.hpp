// SPDX-License-Identifier: MIT

#pragma once

#include <utility>
#include <vector>

#include <geometry.hpp>
#include <irighthandside.hpp>

#include "iboltzmannsolver.hpp"

class SplitRightHandSideSolver : public IBoltzmannSolver
{
    IBoltzmannSolver const& m_boltzmann_solver;

    std::vector<std::reference_wrapper<IRightHandSide const>> m_rhs;

public:
    SplitRightHandSideSolver(
            IBoltzmannSolver const& vlasov_solver,
            std::vector<std::reference_wrapper<IRightHandSide const>> rhs);

    ~SplitRightHandSideSolver() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, DViewX electric_field, double dt) const override;
};
