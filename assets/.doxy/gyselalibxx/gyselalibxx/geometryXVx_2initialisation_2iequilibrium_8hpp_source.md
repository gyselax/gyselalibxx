

# File iequilibrium.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**iequilibrium.hpp**](geometryXVx_2initialisation_2iequilibrium_8hpp.md)

[Go to the documentation of this file](geometryXVx_2initialisation_2iequilibrium_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    virtual DFieldSpVx operator()(DFieldSpVx allfequilibrium) const = 0;
};
```


