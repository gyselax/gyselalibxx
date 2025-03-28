

# File nulladvectionvx.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**nulladvectionvx.hpp**](nulladvectionvx_8hpp.md)

[Go to the documentation of this file](nulladvectionvx_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"
#include "iadvectionvx.hpp"

template <class IdxRangeFdistribu, class IdxRangeSpatial>
class NullAdvectionVelocity : public IAdvectionV<IdxRangeFdistribu, IdxRangeSpatial>
{
public:
    NullAdvectionVelocity() = default;

    ~NullAdvectionVelocity() override = default;

    DField<IdxRangeFdistribu> operator()(
            DField<IdxRangeFdistribu> allfdistribu,
            [[maybe_unused]] DConstField<IdxRangeSpatial> electric_field,
            [[maybe_unused]] double dt) const override
    {
        return allfdistribu;
    }
};
```


