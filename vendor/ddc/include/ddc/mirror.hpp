// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>

#include "ddc/chunk_span.hpp"
#include "ddc/kokkos_allocator.hpp"

namespace ddc {

/// Returns a new `Chunk` with the same layout as `src` allocated on the memory space `Space::memory_space`.
template <class Space, class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror(
        [[maybe_unused]] Space const& space,
        ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    return Chunk<
            std::remove_const_t<ElementType>,
            Support,
            KokkosAllocator<std::remove_const_t<ElementType>, typename Space::memory_space>>(
            src.domain());
}

/// Returns a new host `Chunk` with the same layout as `src`.
template <class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror(ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    return create_mirror(Kokkos::DefaultHostExecutionSpace(), src);
}

/// Returns a new `Chunk` with the same layout as `src` allocated on the memory space `Space::memory_space` and operates a deep copy between the two.
template <class Space, class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_and_copy(
        [[maybe_unused]] Space const& space,
        ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    Chunk<std::remove_const_t<ElementType>,
          Support,
          KokkosAllocator<std::remove_const_t<ElementType>, typename Space::memory_space>>
            chunk(src.domain());
    deepcopy(chunk, src);
    return chunk;
}

/// Returns a new host `Chunk` with the same layout as `src` and operates a deep copy between the two.
template <class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_and_copy(ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    return create_mirror_and_copy(Kokkos::DefaultHostExecutionSpace(), src);
}

/// If `src` is accessible from `space` then returns a copy of `src`,
/// otherwise returns a new `Chunk` with the same layout as `src` allocated on the memory space `Space::memory_space`.
template <class Space, class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_view(
        Space const& space,
        ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    if constexpr (Kokkos::SpaceAccessibility<Space, MemorySpace>::accessible) {
        return src;
    } else {
        return create_mirror(space, src);
    }
}

/// If `src` is host accessible then returns a copy of `src`,
/// otherwise returns a new host `Chunk` with the same layout.
template <class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_view(ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    return create_mirror_view(Kokkos::DefaultHostExecutionSpace(), src);
}

/// If `src` is accessible from `space` then returns a copy of `src`,
/// otherwise returns a new `Chunk` with the same layout as `src` allocated on the memory space `Space::memory_space` and operates a deep copy between the two.
template <class Space, class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_view_and_copy(
        Space const& space,
        ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    if constexpr (Kokkos::SpaceAccessibility<Space, MemorySpace>::accessible) {
        return src;
    } else {
        return create_mirror_and_copy(space, src);
    }
}

/// If `src` is host accessible then returns a copy of `src`,
/// otherwise returns a new host `Chunk` with the same layout as `src` and operates a deep copy between the two.
template <class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_view_and_copy(ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    return create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), src);
}

} // namespace ddc
