// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>

#include "ddc/chunk_span.hpp"

namespace ddc {

/** Copy the content of a borrowed chunk into another
 * @param[out] dst the borrowed chunk in which to copy
 * @param[in]  src the borrowed chunk from which to copy
 * @return dst as a ChunkSpan
*/
template <class ChunkDst, class ChunkSrc>
auto deepcopy(ChunkDst&& dst, ChunkSrc&& src)
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

} // namespace ddc
