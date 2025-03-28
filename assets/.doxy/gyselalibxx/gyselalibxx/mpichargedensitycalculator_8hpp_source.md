

# File mpichargedensitycalculator.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**poisson**](dir_14c5eb4d397dfd4e1a4d5c7bede9e118.md) **>** [**mpichargedensitycalculator.hpp**](mpichargedensitycalculator_8hpp.md)

[Go to the documentation of this file](mpichargedensitycalculator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <mpi.h>

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "geometry.hpp"
#include "ichargedensitycalculator.hpp"
#include "quadrature.hpp"

class MpiChargeDensityCalculator : public IChargeDensityCalculator
{
private:
    IChargeDensityCalculator const& m_local_charge_density_calculator;
    MPI_Comm m_comm;

public:
    explicit MpiChargeDensityCalculator(
            MPI_Comm comm,
            IChargeDensityCalculator const& local_charge_density_calculator);

    void operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const final;
};
```


