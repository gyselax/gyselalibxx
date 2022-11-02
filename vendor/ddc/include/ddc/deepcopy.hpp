// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>
#include <utility>

#include <Kokkos_Core.hpp>

#include "ddc/chunk_span.hpp"
#include "ddc/for_each.hpp"

namespace ddc {

/** Copy the content of a view into another
 * @param[out] dst the view in which to copy
 * @param[in]  src the view from which to copy
 * @return to
 */
template <class ChunkDst, class ChunkSrc>
inline ChunkDst const& deepcopy(ChunkDst&& dst, ChunkSrc&& src) noexcept
{
    static_assert(is_chunk_v<ChunkDst>);
    static_assert(is_chunk_v<ChunkSrc>);
    static_assert(
            std::is_assignable_v<decltype(*dst.data()), decltype(*src.data())>,
            "Not assignable");
    assert(dst.domain().extents() == src.domain().extents());
    Kokkos::deep_copy(dst.allocation_kokkos_view(), src.allocation_kokkos_view());
    return dst;
}

} // namespace ddc
