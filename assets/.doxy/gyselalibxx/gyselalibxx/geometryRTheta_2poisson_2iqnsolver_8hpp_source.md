

# File iqnsolver.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**poisson**](dir_131fdd0509f46f459997bddabd4481b1.md) **>** [**iqnsolver.hpp**](geometryRTheta_2poisson_2iqnsolver_8hpp.md)

[Go to the documentation of this file](geometryRTheta_2poisson_2iqnsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"

class IQNSolver
{
public:
    virtual ~IQNSolver() = default;

    virtual void operator()(
            host_t<DFieldRTheta> electrostatic_potential,
            host_t<DVectorFieldRTheta<X, Y>> electric_field,
            host_t<DConstFieldRTheta> allfdistribu) const = 0;
};
```


