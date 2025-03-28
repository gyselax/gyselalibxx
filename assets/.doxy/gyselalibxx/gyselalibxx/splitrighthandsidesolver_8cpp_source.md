

# File splitrighthandsidesolver.cpp

[**File List**](files.md) **>** [**boltzmann**](dir_7559acab695a99e26dbd57f46ed1b0cd.md) **>** [**splitrighthandsidesolver.cpp**](splitrighthandsidesolver_8cpp.md)

[Go to the documentation of this file](splitrighthandsidesolver_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include <utility>
#include <vector>

#include "geometry.hpp"
#include "iboltzmannsolver.hpp"
#include "irighthandside.hpp"
#include "splitrighthandsidesolver.hpp"

SplitRightHandSideSolver::SplitRightHandSideSolver(
        IBoltzmannSolver const& boltzmann_solver,
        std::vector<std::reference_wrapper<IRightHandSide const>> rhs)
    : m_boltzmann_solver(boltzmann_solver)
    , m_rhs(std::move(rhs))
{
}

DFieldSpXVx SplitRightHandSideSolver::operator()(
        DFieldSpXVx const allfdistribu,
        DConstFieldX const electric_field,
        double const dt) const
{
    for (auto rhsit = m_rhs.begin(); rhsit != m_rhs.end(); ++rhsit) {
        (*rhsit)(allfdistribu, dt / 2.);
    }
    m_boltzmann_solver(allfdistribu, electric_field, dt);
    for (auto rhsit = m_rhs.rbegin(); rhsit != m_rhs.rend(); ++rhsit) {
        (*rhsit)(allfdistribu, dt / 2.);
    }

    return allfdistribu;
}
```


