

# File ddc\_aliases.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**ddc\_aliases.hpp**](ddc__aliases_8hpp.md)

[Go to the documentation of this file](ddc__aliases_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

template <class... Dims>
using Coord = ddc::Coordinate<Dims...>;

template <class... GridTypes>
using Idx = ddc::DiscreteElement<GridTypes...>;

template <class... GridTypes>
using IdxStep = ddc::DiscreteVector<GridTypes...>;

template <class... GridTypes>
using IdxRange = ddc::DiscreteDomain<GridTypes...>;

template <
        class ElementType,
        class IdxRange,
        class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
using FieldMem = ddc::Chunk<ElementType, IdxRange, ddc::KokkosAllocator<ElementType, MemSpace>>;

template <class IdxRange, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
using DFieldMem = FieldMem<double, IdxRange, MemSpace>;

template <
        class ElementType,
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using Field = ddc::ChunkSpan<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;

template <
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DField = Field<double, IdxRange, MemorySpace, LayoutStridedPolicy>;

template <
        class ElementType,
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using ConstField = ddc::ChunkView<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;

template <
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DConstField = ConstField<double, IdxRange, MemorySpace, LayoutStridedPolicy>;

template <class Dim>
using UniformGridBase = ddc::UniformPointSampling<Dim>;

template <class Dim>
using NonUniformGridBase = ddc::NonUniformPointSampling<Dim>;
```


