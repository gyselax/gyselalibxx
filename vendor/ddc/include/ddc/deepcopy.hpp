#pragma once

#include "ddc/block_span.hpp"

namespace detail {
template <class ElementType, class... Meshes, class Layout, class Functor, class... MCoords>
inline void for_each_impl(
        const BlockSpan<ElementType, ProductMDomain<Meshes...>, Layout>& to,
        Functor&& f,
        MCoords&&... mcoords) noexcept
{
    if constexpr (
            sizeof...(MCoords)
            == BlockSpan<ElementType, ProductMDomain<Meshes...>, Layout>::rank()) {
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
        class ElementType,
        class OElementType,
        class... Meshes,
        class... OMeshes,
        class Layout,
        class OLayout>
inline BlockSpan<ElementType, ProductMDomain<Meshes...>, Layout> const& deepcopy(
        BlockSpan<ElementType, ProductMDomain<Meshes...>, Layout> const& to,
        BlockSpan<OElementType, ProductMDomain<OMeshes...>, OLayout> const& from) noexcept
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
template <class ElementType, class... Meshes, class Layout, class Functor>
inline void for_each(
        const BlockSpan<ElementType, ProductMDomain<Meshes...>, Layout>& view,
        Functor&& f) noexcept
{
    detail::for_each_impl(view, std::forward<Functor>(f));
}
