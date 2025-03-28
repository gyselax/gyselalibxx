

# File iinitialisation.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**iinitialisation.hpp**](geometryXVx_2initialisation_2iinitialisation_8hpp.md)

[Go to the documentation of this file](geometryXVx_2initialisation_2iinitialisation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class IInitialisation
{
public:
    virtual ~IInitialisation() = default;

    virtual DFieldSpXVx operator()(DFieldSpXVx allfdistribu) const = 0;
};
```


