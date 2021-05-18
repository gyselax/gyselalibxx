#pragma once

#include <view.h>

#include "mcoord.h"
#include "rcoord.h"

template <class... Tags>
class NonUniformMesh;

namespace detail {

template <
        class... SelectedTags,
        class TagsHead,
        class... TagsQueue,
        class... AllTags,
        class SliceSpec>
inline constexpr auto append_if_all(
        NonUniformMesh<SelectedTags...>,
        NonUniformMesh<TagsHead>,
        NonUniformMesh<AllTags...> m,
        SliceSpec) noexcept
{
    static_assert(
            std::is_integral_v<
                    SliceSpec> || std::is_same_v<std::experimental::all_type, SliceSpec>);
    if constexpr (std::is_integral_v<SliceSpec>) {
        return NonUniformMesh<SelectedTags...>(m);
    } else {
        return NonUniformMesh<TagsHead, SelectedTags...>(m);
    }
}


template <class TagsHead, class... TagsQueue, class SliceSpecsHead, class... SliceSpecsQueue>
inline constexpr auto select_tags(
        NonUniformMesh<TagsHead, TagsQueue...> m,
        SliceSpecsHead&& h,
        SliceSpecsQueue&&... q) noexcept
{
    static_assert(sizeof...(TagsQueue) == sizeof...(SliceSpecsQueue));
    if constexpr (sizeof...(TagsQueue) > 0) {
        return append_if_all(
                select_tags(NonUniformMesh<TagsQueue...>(m), q...),
                NonUniformMesh<TagsHead>(m),
                m,
                h);
    } else {
        return append_if_all(NonUniformMesh<>(RCoord<>(), RLength<>()), m, m, h);
    }
}

} // namespace detail

template <class Tag>
class NonUniformMesh
{
public:
    using RCoord_ = RCoord<Tag>;

    using RLength_ = RLength<Tag>;

    using MCoord_ = MCoord<Tag>;

private:
    /// origin
    TaggedArray<RCoord_> m_origin;

    /// step size
    RLength_ m_step;

    template <class>
    friend class NonUniformMesh;

public:
    template <class... OTags>
    inline constexpr NonUniformMesh(const NonUniformMesh<OTags...>& other) noexcept
        : m_origin(other.m_origin)
        , m_step(other.m_step)
    {
    }

    template <class... OTags>
    inline constexpr NonUniformMesh(NonUniformMesh<OTags...>&& other) noexcept
        : m_origin(std::move(other.m_origin))
        , m_step(std::move(other.m_step))
    {
    }

    template <class OriginType, class StepType>
    inline constexpr NonUniformMesh(OriginType&& origin, StepType&& step) noexcept
        : m_origin(std::forward<OriginType>(origin))
        , m_step(std::forward<StepType>(step))
    {
    }

    friend constexpr bool operator==(const NonUniformMesh& xx, const NonUniformMesh& yy)
    {
        return (&xx == &yy) || (xx.m_origin == yy.m_origin && xx.m_step == yy.m_step);
    }

    friend constexpr bool operator!=(const NonUniformMesh& xx, const NonUniformMesh& yy)
    {
        return (&xx != &yy) && (xx.m_origin != yy.m_origin || xx.m_step != yy.m_step);
    }

    template <class... OTags>
    friend constexpr bool operator==(const NonUniformMesh& xx, const NonUniformMesh<OTags...>& yy)
    {
        return false;
    }

    template <class... OTags>
    friend constexpr bool operator!=(const NonUniformMesh& xx, const NonUniformMesh<OTags...>& yy)
    {
        return false;
    }

    static inline constexpr size_t rank() noexcept
    {
        return sizeof...(Tags);
    }

    inline constexpr RCoord_ origin() const noexcept
    {
        return m_origin;
    }

    inline constexpr RLength_ step() const noexcept
    {
        return m_step;
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> to_real(const MCoord<OTags...> icoord) const noexcept
    {
        return RCoord_((::get<OTags>(origin()) + ::get<OTags>(icoord) * ::get<OTags>(m_step))...);
    }
};

template <class... Tags>
std::ostream& operator<<(std::ostream& out, NonUniformMesh<Tags...> const& dom)
{
    return out << "NonUniformMesh( origin=" << dom.origin() << ", unitvec=" << dom.step() << " )";
}

template <class... Tags, class... SliceSpecs>
inline constexpr auto submesh(const NonUniformMesh<Tags...>& mesh, SliceSpecs... slices) noexcept
{
    using ReturnType = decltype(detail::select_tags(mesh, std::forward<SliceSpecs>(slices)...));
    return ReturnType(mesh);
}

template <class... Tags>
using Mesh = NonUniformMesh<Tags...>;

using MeshX = Mesh<Dim::X>;

using MeshVx = Mesh<Dim::Vx>;

using MeshXVx = Mesh<Dim::X, Dim::Vx>;
