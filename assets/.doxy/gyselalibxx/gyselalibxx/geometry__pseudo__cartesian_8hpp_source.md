

# File geometry\_pseudo\_cartesian.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**geometry\_pseudo\_cartesian.hpp**](geometry__pseudo__cartesian_8hpp.md)

[Go to the documentation of this file](geometry__pseudo__cartesian_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


struct X_pC
{
    static bool constexpr PERIODIC = false;

    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = X_pC;
};

struct Y_pC
{
    static bool constexpr PERIODIC = false;

    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = Y_pC;
};

using CoordX_pC = Coord<X_pC>;
using CoordY_pC = Coord<Y_pC>;
using CoordXY_pC = Coord<X_pC, Y_pC>;
```


