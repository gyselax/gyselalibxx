

# File predcorr.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**time\_integration**](dir_61df5a05e7e6762880c4f92d6d795362.md) **>** [**predcorr.hpp**](geometryXVx_2time__integration_2predcorr_8hpp.md)

[Go to the documentation of this file](geometryXVx_2time__integration_2predcorr_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "itimesolver.hpp"

class IQNSolver;
class IBoltzmannSolver;

class PredCorr : public ITimeSolver
{
private:
    IBoltzmannSolver const& m_boltzmann_solver;

    IQNSolver const& m_poisson_solver;

public:
    PredCorr(IBoltzmannSolver const& boltzmann_solver, IQNSolver const& poisson_solver);

    ~PredCorr() override = default;

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double time_start, double dt, int steps = 1)
            const override;
};
```


