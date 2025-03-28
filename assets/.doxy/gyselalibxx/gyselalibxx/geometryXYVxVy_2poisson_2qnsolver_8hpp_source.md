

# File qnsolver.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**poisson**](dir_14c5eb4d397dfd4e1a4d5c7bede9e118.md) **>** [**qnsolver.hpp**](geometryXYVxVy_2poisson_2qnsolver_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2poisson_2qnsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include "chargedensitycalculator.hpp"
#include "ddc_aliases.hpp"
#include "ipoisson_solver.hpp"
#include "iqnsolver.hpp"

class QNSolver : public IQNSolver
{
    using PoissonSolver = IPoissonSolver<
            IdxRangeXY,
            IdxRangeXY,
            typename Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::layout_right>;
    PoissonSolver const& m_solve_poisson;
    IChargeDensityCalculator const& m_compute_rho;

public:
    QNSolver(PoissonSolver const& solve_poisson, IChargeDensityCalculator const& compute_rho);

    ~QNSolver() override = default;

    void operator()(
            DFieldXY electrostatic_potential,
            DFieldXY electric_field_x,
            DFieldXY electric_field_y,
            DConstFieldSpVxVyXY allfdistribu) const override;
};
```


