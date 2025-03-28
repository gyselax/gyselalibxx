

# File iadvectionvx.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**iadvectionvx.hpp**](iadvectionvx_8hpp.md)

[Go to the documentation of this file](iadvectionvx_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

template <class Geometry, class GridV>
class IAdvectionVelocity
{
public:
    virtual ~IAdvectionVelocity() = default;

    virtual DField<typename Geometry::IdxRangeFdistribu> operator()(
            DField<typename Geometry::IdxRangeFdistribu> allfdistribu,
            DConstField<typename Geometry::IdxRangeSpatial> electric_field,
            double dt) const = 0;
};
```


