// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>
#include <utility>

#include "ddc/chunk_span.hpp"
#include "ddc/for_each.hpp"

/** Copy the content of a view into another
 * @param[out] to    the view in which to copy
 * @param[in]  from  the view from which to copy
 * @return to
 */
template <class ChunkDst, class ChunkSrc>
inline ChunkDst const& deepcopy(ChunkDst&& to, ChunkSrc&& from) noexcept
{
    static_assert(is_chunk_v<ChunkDst>);
    static_assert(is_chunk_v<ChunkSrc>);
    static_assert(
            std::is_assignable_v<decltype(*to.data()), decltype(*from.data())>,
            "Not assignable");
    assert(to.domain().extents() == from.domain().extents());
    for_each_n(to.domain().extents(), [&to, &from](auto&& idx) {
        to(to.domain().front() + idx) = from(from.domain().front() + idx);
    });
    return to;
}
