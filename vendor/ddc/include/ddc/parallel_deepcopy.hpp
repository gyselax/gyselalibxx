// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>

#include "ddc/chunk_common.hpp"

namespace ddc {

/** Copy the content of a borrowed chunk into another
 * @param[out] dst the borrowed chunk in which to copy
 * @param[in]  src the borrowed chunk from which to copy
 * @return dst as a ChunkSpan
*/
template <class ChunkDst, class ChunkSrc>
auto parallel_deepcopy(ChunkDst&& dst, ChunkSrc&& src)
{
    static_assert(is_borrowed_chunk_v<ChunkDst>);
    static_assert(is_borrowed_chunk_v<ChunkSrc>);
    static_assert(
            std::is_assignable_v<chunk_reference_t<ChunkDst>, chunk_reference_t<ChunkSrc>>,
            "Not assignable");
    assert(dst.domain().extents() == src.domain().extents());
    Kokkos::deep_copy(dst.allocation_kokkos_view(), src.allocation_kokkos_view());
    return dst.span_view();
}

/** Copy the content of a borrowed chunk into another
 * @param[in] execution_space a Kokkos execution space where the loop will be executed on
 * @param[out] dst the borrowed chunk in which to copy
 * @param[in]  src the borrowed chunk from which to copy
 * @return dst as a ChunkSpan
*/
template <class ExecSpace, class ChunkDst, class ChunkSrc>
auto parallel_deepcopy(ExecSpace const& execution_space, ChunkDst&& dst, ChunkSrc&& src)
{
    static_assert(is_borrowed_chunk_v<ChunkDst>);
    static_assert(is_borrowed_chunk_v<ChunkSrc>);
    static_assert(
            std::is_assignable_v<chunk_reference_t<ChunkDst>, chunk_reference_t<ChunkSrc>>,
            "Not assignable");
    assert(dst.domain().extents() == src.domain().extents());
    Kokkos::deep_copy(execution_space, dst.allocation_kokkos_view(), src.allocation_kokkos_view());
    return dst.span_view();
}

} // namespace ddc
