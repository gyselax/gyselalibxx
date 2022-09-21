// SPDX-License-Identifier: MIT

#include <utility>
#include <vector>

#include <geometry.hpp>
#include <isource.hpp>

#include "iboltzmannsolver.hpp"
#include "splitsourceenvironmentsolver.hpp"

SplitSourceEnvironmentSolver::SplitSourceEnvironmentSolver(
        IBoltzmannSolver const& boltzmann_solver,
        std::vector<std::reference_wrapper<ISource const>> rhs)
    : m_boltzmann_solver(boltzmann_solver)
    , m_rhs(std::move(rhs))
{
}

DSpanSpXVx SplitSourceEnvironmentSolver::operator()(
        DSpanSpXVx const allfdistribu,
        DViewX const electric_field,
        double const dt) const
{
    for (auto rhsit = m_rhs.begin(); rhsit != m_rhs.end(); ++rhsit) {
        (*rhsit)(allfdistribu, dt);
    }
    m_boltzmann_solver(allfdistribu, electric_field, dt);
    for (auto rhsit = m_rhs.rbegin(); rhsit != m_rhs.rend(); ++rhsit) {
        (*rhsit)(allfdistribu, dt);
    }

    return allfdistribu;
}
