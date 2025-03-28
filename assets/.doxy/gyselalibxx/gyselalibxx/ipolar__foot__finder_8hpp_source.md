

# File ipolar\_foot\_finder.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**ipolar\_foot\_finder.hpp**](ipolar__foot__finder_8hpp.md)

[Go to the documentation of this file](ipolar__foot__finder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include "ddc_aliases.hpp"
#include "vector_field.hpp"

template <
        class GridRadial,
        class GridPoloidal,
        class AdvectionDim1,
        class AdvectionDim2,
        class MemorySpace>
class IPolarFootFinder
{
protected:
    using GridR = GridRadial;
    using GridTheta = GridPoloidal;

    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;

    using X = AdvectionDim1;
    using Y = AdvectionDim2;

    using memory_space = MemorySpace;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;

public:
    virtual ~IPolarFootFinder() = default;

    virtual void operator()(
            Field<Coord<R, Theta>, IdxRangeRTheta, memory_space> feet,
            DVectorConstField<IdxRangeRTheta, VectorIndexSet<X, Y>, memory_space> advection_field,
            double dt) const = 0;
};
```


