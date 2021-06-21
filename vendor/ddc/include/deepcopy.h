#pragma once

#include "blockview.h"

namespace detail {
template <class Mesh, class ElementType, bool CONTIGUOUS, class Functor, class... Indices>
inline void for_each_impl(
        const BlockView<MDomainImpl<Mesh>, ElementType, CONTIGUOUS>& to,
        Functor&& f,
        Indices&&... idxs) noexcept
{
    if constexpr (
            sizeof...(Indices) == BlockView<MDomainImpl<Mesh>, ElementType, CONTIGUOUS>::rank()) {
        f(std::forward<Indices>(idxs)...);
    } else {
        for (std::size_t ii = 0; ii < to.extent(sizeof...(Indices)); ++ii) {
            for_each_impl(to, std::forward<Functor>(f), std::forward<Indices>(idxs)..., ii);
        }
    }
}
template <class MCoord>
struct sequential_for_impl
{
    template <class Extents, class Functor, class... Indices>
    inline void operator()(const Extents& extents, Functor&& f, Indices&&... idxs) const noexcept
    {
        if constexpr (sizeof...(Indices) == Extents::rank()) {
            f(MCoord(std::forward<Indices>(idxs)...));
        } // namespace detail
        else {
            for (std::size_t ii = 0; ii < extents.extent(sizeof...(Indices)); ++ii) {
                (*this)(extents, std::forward<Functor>(f), std::forward<Indices>(idxs)..., ii);
            }
        }
    }
};
} // namespace detail

/** Copy the content of a view into another
 * @param[out] to    the view in which to copy
 * @param[in]  from  the view from which to copy
 * @return to
 */
template <class... Tags, class ElementType, bool CONTIGUOUS, bool OCONTIGUOUS>
inline BlockView<UniformMDomain<Tags...>, ElementType, CONTIGUOUS> deepcopy(
        BlockView<UniformMDomain<Tags...>, ElementType, CONTIGUOUS> to,
        BlockView<UniformMDomain<Tags...>, ElementType, OCONTIGUOUS> const& from) noexcept
{
    assert(to.extents() == from.extents());
    for_each(to, [&to, &from](auto... idxs) { to(idxs...) = from(idxs...); });
    return to;
}

/** Copy the content of a view into another
 * @param[out] to    the view in which to copy
 * @param[in]  from  the view from which to copy
 * @return to
 */
template <
        class... Tags,
        class... OTags,
        class ElementType,
        class OElementType,
        bool CONTIGUOUS,
        bool OCONTIGUOUS>
inline BlockView<UniformMDomain<Tags...>, ElementType, CONTIGUOUS> deepcopy(
        BlockView<UniformMDomain<Tags...>, ElementType, CONTIGUOUS> to,
        BlockView<UniformMDomain<OTags...>, OElementType, OCONTIGUOUS> const& from) noexcept
{
    static_assert(std::is_convertible_v<OElementType, ElementType>, "Not convertible");
    using MCoord_ = typename UniformMDomain<Tags...>::mcoord_type;
    assert(to.extents() == from.extents());
    constexpr auto sequential_for = detail::sequential_for_impl<MCoord_>();
    sequential_for(to.extents(), [&to, &from](const auto& domain) { to(domain) = from(domain); });
    return to;
}

/** iterates over the domain of a view
 * @param[in] view  the view whose domain to iterate
 * @param[in] f     a functor taking the list of indices as parameter
 */
template <class... Tags, class ElementType, bool CONTIGUOUS, class Functor>
inline void for_each(
        const BlockView<UniformMDomain<Tags...>, ElementType, CONTIGUOUS>& view,
        Functor&& f) noexcept
{
    detail::for_each_impl(view, std::forward<Functor>(f));
}
