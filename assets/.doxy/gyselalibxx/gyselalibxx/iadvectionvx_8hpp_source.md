

# File iadvectionvx.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**iadvectionvx.hpp**](iadvectionvx_8hpp.md)

[Go to the documentation of this file](iadvectionvx_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

template <class Geometry, class GridV, class DataType = double>
class IAdvectionVelocity
{
public:
    virtual ~IAdvectionVelocity() = default;

    virtual Field<DataType, typename Geometry::IdxRangeFdistribu> operator()(
            Field<DataType, typename Geometry::IdxRangeFdistribu> allfdistribu,
            ConstField<DataType, typename Geometry::IdxRangeSpatial> electric_field,
            DataType dt) const = 0;
};
```


