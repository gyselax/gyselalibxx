

# File ichargedensitycalculator.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**poisson**](dir_d78fdb6d05340e24a2e187de33ea09a4.md) **>** [**ichargedensitycalculator.hpp**](geometryXVx_2poisson_2ichargedensitycalculator_8hpp.md)

[Go to the documentation of this file](geometryXVx_2poisson_2ichargedensitycalculator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"

class IChargeDensityCalculator
{
public:
    virtual DFieldX operator()(DFieldX rho, DConstFieldSpXVx allfdistribu) const = 0;
};
```


