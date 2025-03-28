

# File iqnsolver.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**poisson**](dir_14c5eb4d397dfd4e1a4d5c7bede9e118.md) **>** [**iqnsolver.hpp**](geometryXYVxVy_2poisson_2iqnsolver_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2poisson_2iqnsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

class IQNSolver
{
public:
    virtual ~IQNSolver() = default;

    virtual void operator()(
            DFieldXY electrostatic_potential,
            DFieldXY electric_field_x,
            DFieldXY electric_field_y,
            DConstFieldSpVxVyXY allfdistribu) const = 0;
};
```


