

# File ichargedensitycalculator.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**poisson**](dir_14c5eb4d397dfd4e1a4d5c7bede9e118.md) **>** [**ichargedensitycalculator.hpp**](geometryXYVxVy_2poisson_2ichargedensitycalculator_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2poisson_2ichargedensitycalculator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"

class IChargeDensityCalculator
{
public:
    virtual void operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const = 0;
};
```


