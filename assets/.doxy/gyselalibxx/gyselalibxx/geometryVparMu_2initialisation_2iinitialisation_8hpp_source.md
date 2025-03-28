

# File iinitialisation.hpp

[**File List**](files.md) **>** [**geometryVparMu**](dir_9a2f28dc8f538ee0f4428810facf29b8.md) **>** [**initialisation**](dir_99d29839093a8e7b0be0d596be7efa54.md) **>** [**iinitialisation.hpp**](geometryVparMu_2initialisation_2iinitialisation_8hpp.md)

[Go to the documentation of this file](geometryVparMu_2initialisation_2iinitialisation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class IInitialisation
{
public:
    virtual ~IInitialisation() = default;

    virtual DFieldSpVparMu operator()(DFieldSpVparMu allfdistribu) const = 0;
};
```


