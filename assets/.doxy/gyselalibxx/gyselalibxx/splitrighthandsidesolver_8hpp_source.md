

# File splitrighthandsidesolver.hpp

[**File List**](files.md) **>** [**boltzmann**](dir_7559acab695a99e26dbd57f46ed1b0cd.md) **>** [**splitrighthandsidesolver.hpp**](splitrighthandsidesolver_8hpp.md)

[Go to the documentation of this file](splitrighthandsidesolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <utility>
#include <vector>

#include "geometry.hpp"
#include "iboltzmannsolver.hpp"
#include "irighthandside.hpp"

class SplitRightHandSideSolver : public IBoltzmannSolver
{
    IBoltzmannSolver const& m_boltzmann_solver;

    std::vector<std::reference_wrapper<IRightHandSide const>> m_rhs;

public:
    SplitRightHandSideSolver(
            IBoltzmannSolver const& vlasov_solver,
            std::vector<std::reference_wrapper<IRightHandSide const>> rhs);

    ~SplitRightHandSideSolver() override = default;

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, DConstFieldX electric_field, double dt)
            const override;
};
```


