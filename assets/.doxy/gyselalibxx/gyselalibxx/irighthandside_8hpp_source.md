

# File irighthandside.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**irighthandside.hpp**](irighthandside_8hpp.md)

[Go to the documentation of this file](irighthandside_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

enum class RhsType { Source, Sink };

class IRightHandSide
{
public:
    virtual ~IRightHandSide() = default;

    virtual DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const = 0;
};
```


