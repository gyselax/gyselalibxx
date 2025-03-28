

# File iboltzmannsolver.hpp

[**File List**](files.md) **>** [**boltzmann**](dir_7559acab695a99e26dbd57f46ed1b0cd.md) **>** [**iboltzmannsolver.hpp**](iboltzmannsolver_8hpp.md)

[Go to the documentation of this file](iboltzmannsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class IBoltzmannSolver
{
public:
    virtual ~IBoltzmannSolver() = default;

    virtual DFieldSpXVx operator()(DFieldSpXVx allfdistribu, DConstFieldX efield, double dt)
            const = 0;
};
```


