

# File edge.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**edge.hpp**](edge_8hpp.md)

[Go to the documentation of this file](edge_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry_descriptors.hpp"


template <class Patch, class Grid1D, Extremity extremity_val>
struct Edge
{
    using associated_patch = Patch;
    using perpendicular_grid = Grid1D;
    using parallel_grid = std::conditional_t<
            std::is_same_v<Grid1D, typename Patch::Grid1>,
            typename Patch::Grid2,
            typename Patch::Grid1>;
    static constexpr Extremity extremity = extremity_val;
};


struct OutsideEdge
{
};
```


