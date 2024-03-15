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
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
    return Chunk(
            src.domain(),
            KokkosAllocator<std::remove_const_t<ElementType>, typename Space::memory_space>());
}

/// Returns a new host `Chunk` with the same layout as `src`.
template <class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror(ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
    return Chunk(src.domain(), HostAllocator<std::remove_const_t<ElementType>>());
}

/// Returns a new `Chunk` with the same layout as `src` allocated on the memory space `Space::memory_space` and operates a deep copy between the two.
template <class Space, class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_and_copy(
        Space const& space,
        ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
    Chunk chunk = create_mirror(space, src);
    parallel_deepcopy(space, chunk, src);
    return chunk;
}

/// Returns a new host `Chunk` with the same layout as `src` and operates a deep copy between the two.
template <class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_and_copy(ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
    Chunk chunk = create_mirror(src);
    parallel_deepcopy(chunk, src);
    return chunk;
}

/// If `src` is accessible from `space` then returns a copy of `src`,
/// otherwise returns a new `Chunk` with the same layout as `src` allocated on the memory space `Space::memory_space`.
template <class Space, class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_view(
        Space const& space,
        ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
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
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::HostSpace, MemorySpace>::accessible) {
        return src;
    } else {
        return create_mirror(src);
    }
}

/// If `src` is accessible from `space` then returns a copy of `src`,
/// otherwise returns a new `Chunk` with the same layout as `src` allocated on the memory space `Space::memory_space` and operates a deep copy between the two.
template <class Space, class ElementType, class Support, class Layout, class MemorySpace>
auto create_mirror_view_and_copy(
        Space const& space,
        ChunkSpan<ElementType, Support, Layout, MemorySpace> const& src)
{
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
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
    static_assert(
            std::is_same_v<Layout, std::experimental::layout_right>,
            "Only layout right is supported");
    if constexpr (Kokkos::SpaceAccessibility<Kokkos::HostSpace, MemorySpace>::accessible) {
        return src;
    } else {
        return create_mirror_and_copy(src);
    }
}

} // namespace ddc
