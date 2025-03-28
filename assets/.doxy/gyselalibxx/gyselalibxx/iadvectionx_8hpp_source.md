

# File iadvectionx.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**iadvectionx.hpp**](iadvectionx_8hpp.md)

[Go to the documentation of this file](iadvectionx_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

template <class Geometry, class GridX>
class IAdvectionSpatial
{
public:
    virtual ~IAdvectionSpatial() = default;
    virtual DField<typename Geometry::IdxRangeFdistribu> operator()(
            DField<typename Geometry::IdxRangeFdistribu> allfdistribu,
            double dt) const = 0;
};
```


