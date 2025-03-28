

# File iadvection\_rtheta.hpp

[**File List**](files.md) **>** [**advection**](dir_18bb63d4be19d3ea733e61d8625caf4d.md) **>** [**iadvection\_rtheta.hpp**](iadvection__rtheta_8hpp.md)

[Go to the documentation of this file](iadvection__rtheta_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

class IAdvectionRTheta
{
public:
    virtual ~IAdvectionRTheta() = default;

    virtual DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<X, Y> advection_field,
            double const dt) const = 0;

    virtual DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<R, Theta> advection_field,
            CoordXY const& advection_field_xy_centre,
            double const dt) const = 0;
};
```


