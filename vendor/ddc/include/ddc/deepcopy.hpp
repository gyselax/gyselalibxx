#pragma once

#include "ddc/block_span.hpp"

namespace detail {
template <class... Meshes, class ElementType, class Layout, class Functor, class... MCoords>
inline void for_each_impl(
        const BlockSpan<ProductMDomain<Meshes...>, ElementType, Layout>& to,
        Functor&& f,
        MCoords&&... mcoords) noexcept
{
    if constexpr (
            sizeof...(MCoords)
            == BlockSpan<ProductMDomain<Meshes...>, ElementType, Layout>::rank()) {
        f(std::forward<MCoords>(mcoords)...);
    } else {
        using CurrentMesh = type_seq_element_t<sizeof...(MCoords), detail::TypeSeq<Meshes...>>;
        for (auto&& ii : get_domain<CurrentMesh>(to)) {
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
template <
        class... Meshes,
        class... OMeshes,
        class ElementType,
        class OElementType,
        class Layout,
        class OLayout>
inline BlockSpan<ProductMDomain<Meshes...>, ElementType, Layout> const& deepcopy(
        BlockSpan<ProductMDomain<Meshes...>, ElementType, Layout> const& to,
        BlockSpan<ProductMDomain<OMeshes...>, OElementType, OLayout> const& from) noexcept
{
    static_assert(std::is_convertible_v<OElementType, ElementType>, "Not convertible");
    assert(to.domain().front() == from.domain().front());
    assert(to.domain().back() == from.domain().back());
    for_each(to, [&to, &from](auto&&... idxs) { to(idxs...) = from(idxs...); });
    return to;
}

/** iterates over the domain of a view
 * @param[in] view  the view whose domain to iterate
 * @param[in] f     a functor taking the list of indices as parameter
 */
template <class... Meshes, class ElementType, class Layout, class Functor>
inline void for_each(
        const BlockSpan<ProductMDomain<Meshes...>, ElementType, Layout>& view,
        Functor&& f) noexcept
{
    detail::for_each_impl(view, std::forward<Functor>(f));
}
