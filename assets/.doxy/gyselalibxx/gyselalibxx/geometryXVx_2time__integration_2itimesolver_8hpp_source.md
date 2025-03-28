

# File itimesolver.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**time\_integration**](dir_61df5a05e7e6762880c4f92d6d795362.md) **>** [**itimesolver.hpp**](geometryXVx_2time__integration_2itimesolver_8hpp.md)

[Go to the documentation of this file](geometryXVx_2time__integration_2itimesolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class ITimeSolver
{
public:
    virtual ~ITimeSolver() = default;

    virtual DFieldSpXVx operator()(
            DFieldSpXVx allfdistribu,
            double time_start,
            double dt,
            int steps = 1) const = 0;
};
```


