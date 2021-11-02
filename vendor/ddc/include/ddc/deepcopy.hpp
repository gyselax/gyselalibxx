#pragma once

#include "ddc/chunk_span.hpp"

namespace detail {
template <class ElementType, class... DDims, class Layout, class Functor, class... MCoords>
inline void for_each_impl(
        ChunkSpan<ElementType, DiscreteDomain<DDims...>, Layout> const& to,
        Functor&& f,
        MCoords&&... mcoords) noexcept
{
    if constexpr (
            sizeof...(MCoords)
            == ChunkSpan<ElementType, DiscreteDomain<DDims...>, Layout>::rank()) {
        f(std::forward<MCoords>(mcoords)...);
    } else {
        using CurrentDDim = type_seq_element_t<sizeof...(MCoords), detail::TypeSeq<DDims...>>;
        for (auto&& ii : get_domain<CurrentDDim>(to)) {
            for_each_impl(to, std::forward<Functor>(f), std::forward<MCoords>(mcoords)..., ii);
        }
    }
}
} // namespace detail

/** Copy the content of a view into another
 * @param[out] to    the view in which to copy
 * @param[in]  from  the view from which to copy
 * @return to
 */
template <class ChunkDst, class ChunkSrc>
inline ChunkDst const& deepcopy(ChunkDst&& to, ChunkSrc&& from) noexcept
{
    static_assert(is_chunkspan_v<ChunkDst>);
    static_assert(is_chunkspan_v<ChunkSrc>);
    static_assert(
            std::is_assignable_v<decltype(*to.data()), decltype(*from.data())>,
            "Not assignable");
    assert(to.domain().front() == from.domain().front());
    assert(to.domain().back() == from.domain().back());
    for_each(to.span_view(), [&to, &from](auto&&... idxs) { to(idxs...) = from(idxs...); });
    return to;
}

/** iterates over the domain of a view
 * @param[in] view  the view whose domain to iterate
 * @param[in] f     a functor taking the list of indices as parameter
 */
template <class ElementType, class... DDims, class Layout, class Functor>
inline void for_each(
        const ChunkSpan<ElementType, DiscreteDomain<DDims...>, Layout>& view,
        Functor&& f) noexcept
{
    detail::for_each_impl(view, std::forward<Functor>(f));
}
