

# File nullqnsolver.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**poisson**](dir_d78fdb6d05340e24a2e187de33ea09a4.md) **>** [**nullqnsolver.hpp**](geometryXVx_2poisson_2nullqnsolver_8hpp.md)

[Go to the documentation of this file](geometryXVx_2poisson_2nullqnsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "iqnsolver.hpp"
class NullQNSolver : public IQNSolver
{
public:
    NullQNSolver() = default;

    ~NullQNSolver() override = default;
    void operator()(
            DFieldX const electrostatic_potential,
            DFieldX const electric_field,
            DConstFieldSpXVx const allfdistribu) const override;
};
```


