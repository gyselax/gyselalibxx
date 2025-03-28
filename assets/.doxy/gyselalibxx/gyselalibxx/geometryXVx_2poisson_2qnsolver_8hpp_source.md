

# File qnsolver.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**poisson**](dir_d78fdb6d05340e24a2e187de33ea09a4.md) **>** [**qnsolver.hpp**](geometryXVx_2poisson_2qnsolver_8hpp.md)

[Go to the documentation of this file](geometryXVx_2poisson_2qnsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"
#include "ichargedensitycalculator.hpp"
#include "ipoisson_solver.hpp"
#include "iqnsolver.hpp"

class QNSolver : public IQNSolver
{
    using PoissonSolver = IPoissonSolver<
            IdxRangeX,
            IdxRangeX,
            typename Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::layout_right>;
    PoissonSolver const& m_solve_poisson;
    IChargeDensityCalculator const& m_compute_rho;

public:
    QNSolver(PoissonSolver const& solve_poisson, IChargeDensityCalculator const& compute_rho);

    ~QNSolver() override = default;

    void operator()(
            DFieldX electrostatic_potential,
            DFieldX electric_field,
            DConstFieldSpXVx allfdistribu) const override;
};
```


