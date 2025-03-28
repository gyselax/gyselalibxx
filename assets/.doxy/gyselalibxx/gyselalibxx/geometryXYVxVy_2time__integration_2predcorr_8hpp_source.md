

# File predcorr.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**time\_integration**](dir_e2479f83d09a2f8b4ff065e45deaef4e.md) **>** [**predcorr.hpp**](geometryXYVxVy_2time__integration_2predcorr_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2time__integration_2predcorr_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "itimesolver.hpp"

class IQNSolver;
class IVlasovSolver;

class PredCorr : public ITimeSolver
{
private:
    IVlasovSolver const& m_vlasov_solver;

    IQNSolver const& m_poisson_solver;

public:
    PredCorr(IVlasovSolver const& vlasov_solver, IQNSolver const& poisson_solver);

    ~PredCorr() override = default;

    DFieldSpVxVyXY operator()(DFieldSpVxVyXY allfdistribu, double dt, int steps = 1) const override;
};
```


