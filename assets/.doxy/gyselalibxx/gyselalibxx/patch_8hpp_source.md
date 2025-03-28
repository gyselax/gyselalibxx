

# File patch.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**patch.hpp**](patch_8hpp.md)

[Go to the documentation of this file](patch_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

template <class... T>
struct Patch;

template <class grid1, class grid2, class bsplines_dim1, class bsplines_dim2>
struct Patch<grid1, grid2, bsplines_dim1, bsplines_dim2>
{
    static int constexpr n_dims = 2;

    using Grid1 = grid1;
    using Grid2 = grid2;

    using Dim1 = typename grid1::continuous_dimension_type;
    using Dim2 = typename grid2::continuous_dimension_type;

    using BSplines1 = bsplines_dim1;
    using BSplines2 = bsplines_dim2;

    using Coord1 = Coord<Dim1>;
    using Coord2 = Coord<Dim2>;
    using Coord12 = Coord<Dim1, Dim2>;

    using Idx1 = Idx<Grid1>;
    using Idx2 = Idx<Grid2>;
    using Idx12 = Idx<Grid1, Grid2>;

    using IdxStep1 = IdxStep<Grid1>;
    using IdxStep2 = IdxStep<Grid2>;
    using IdxStep12 = IdxStep<Grid1, Grid2>;

    using IdxRange1 = IdxRange<Grid1>;
    using IdxRange2 = IdxRange<Grid2>;
    using IdxRange12 = IdxRange<Grid1, Grid2>;

    using IdxRangeBS1 = IdxRange<BSplines1>;
    using IdxRangeBS2 = IdxRange<BSplines2>;
    using IdxRangeBS12 = IdxRange<BSplines1, BSplines2>;
};
```


