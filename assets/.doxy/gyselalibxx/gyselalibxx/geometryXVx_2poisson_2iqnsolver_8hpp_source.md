

# File iqnsolver.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**poisson**](dir_d78fdb6d05340e24a2e187de33ea09a4.md) **>** [**iqnsolver.hpp**](geometryXVx_2poisson_2iqnsolver_8hpp.md)

[Go to the documentation of this file](geometryXVx_2poisson_2iqnsolver_8hpp.md)


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
            DFieldX electrostatic_potential,
            DFieldX electric_field,
            DConstFieldSpXVx allfdistribu) const = 0;
};
```


