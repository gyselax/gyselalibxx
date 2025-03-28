

# File iequilibrium.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**initialisation**](dir_51031f497920158ed20948cdaeaff0bc.md) **>** [**iequilibrium.hpp**](geometryXYVxVy_2initialisation_2iequilibrium_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2initialisation_2iequilibrium_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    virtual DFieldSpVxVy operator()(DFieldSpVxVy allfequilibrium) const = 0;
};
```


