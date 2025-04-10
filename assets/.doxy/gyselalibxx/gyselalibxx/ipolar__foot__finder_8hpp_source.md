

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
        class VectorIndexSetAdvDims,
        class IdxRangeBatched,
        class MemorySpace>
class IPolarFootFinder
{
    // Check types
    static_assert(
            (ddc::is_non_uniform_point_sampling_v<GridRadial>)
            || (ddc::is_uniform_point_sampling_v<GridRadial>));
    static_assert(
            (ddc::is_non_uniform_point_sampling_v<GridPoloidal>)
            || (ddc::is_uniform_point_sampling_v<GridPoloidal>));
    static_assert(is_vector_index_set_v<VectorIndexSetAdvDims>);
    static_assert(ddc::is_discrete_domain_v<IdxRangeBatched>);
    static_assert(Kokkos::is_memory_space_v<MemorySpace>);

    // Check that grids make sense
    static_assert(
            ddc::in_tags_v<GridRadial, ddc::to_type_seq_t<IdxRangeBatched>>,
            "The radial grid must be found in the batched index range");
    static_assert(
            ddc::in_tags_v<GridPoloidal, ddc::to_type_seq_t<IdxRangeBatched>>,
            "The poloidal grid must be found in the batched index range");

    // Check that VectorIndexSetAdvDims makes sense
    static_assert(ddc::type_seq_size_v<VectorIndexSetAdvDims> == 2);

protected:
    using GridR = GridRadial;
    using GridTheta = GridPoloidal;

    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;

    using VectorIndexSetAdvectionDims = VectorIndexSetAdvDims;

public:
    using memory_space = MemorySpace;

    using IdxRangeOperator = IdxRangeBatched;

public:
    virtual ~IPolarFootFinder() = default;

    virtual void operator()(
            Field<Coord<R, Theta>, IdxRangeOperator, memory_space> feet,
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field,
            double dt) const = 0;
};
```


