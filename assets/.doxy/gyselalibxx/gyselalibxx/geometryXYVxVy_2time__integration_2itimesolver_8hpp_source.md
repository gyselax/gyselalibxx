

# File itimesolver.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**time\_integration**](dir_e2479f83d09a2f8b4ff065e45deaef4e.md) **>** [**itimesolver.hpp**](geometryXYVxVy_2time__integration_2itimesolver_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2time__integration_2itimesolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    virtual DFieldSpVxVyXY operator()(DFieldSpVxVyXY allfdistribu, double dt, int steps = 1)
            const = 0;
};
```


