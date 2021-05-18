#pragma once

#include "mcoord.h"
#include "rcoord.h"

template <class... Tags>
class UniformMesh
{
public:
    using RCoord_ = RCoord<Tags...>;

    using RLength_ = RLength<Tags...>;

    using MCoord_ = MCoord<Tags...>;

private:
    /// origin
    RCoord_ m_origin;

    /// step size
    RLength_ m_step;

    template <class...>
    friend class UniformMesh;

public:
    template <class... OTags>
    inline constexpr UniformMesh(const UniformMesh<OTags...>& other) noexcept
        : m_origin(other.m_origin)
        , m_step(other.m_step)
    {
    }

    template <class... OTags>
    inline constexpr UniformMesh(UniformMesh<OTags...>&& other) noexcept
        : m_origin(std::move(other.m_origin))
        , m_step(std::move(other.m_step))
    {
    }

    template <class OriginType, class StepType>
    inline constexpr UniformMesh(OriginType&& origin, StepType&& step) noexcept
        : m_origin(std::forward<OriginType>(origin))
        , m_step(std::forward<StepType>(step))
    {
    }

    friend constexpr bool operator==(const UniformMesh& xx, const UniformMesh& yy)
    {
        return (&xx == &yy) || (xx.m_origin == yy.m_origin && xx.m_step == yy.m_step);
    }

    friend constexpr bool operator!=(const UniformMesh& xx, const UniformMesh& yy)
    {
        return (&xx != &yy) && (xx.m_origin != yy.m_origin || xx.m_step != yy.m_step);
    }

    template <class... OTags>
    friend constexpr bool operator==(const UniformMesh& xx, const UniformMesh<OTags...>& yy)
    {
        return false;
    }

    template <class... OTags>
    friend constexpr bool operator!=(const UniformMesh& xx, const UniformMesh<OTags...>& yy)
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
        return RCoord<OTags...>(
                (::get<OTags>(origin()) + ::get<OTags>(icoord) * ::get<OTags>(m_step))...);
    }
};

namespace detail {

template <
        class... SelectedTags,
        class TagsHead,
        class... TagsQueue,
        class... AllTags,
        class SliceSpec>
inline constexpr auto append_if_all(
        UniformMesh<SelectedTags...>,
        UniformMesh<TagsHead>,
        UniformMesh<AllTags...> m,
        SliceSpec) noexcept
{
    static_assert(
            std::is_integral_v<
                    SliceSpec> || std::is_same_v<std::experimental::all_type, SliceSpec>);
    if constexpr (std::is_integral_v<SliceSpec>) {
        return UniformMesh<SelectedTags...>(m);
    } else {
        return UniformMesh<TagsHead, SelectedTags...>(m);
    }
}


template <class TagsHead, class... TagsQueue, class SliceSpecsHead, class... SliceSpecsQueue>
inline constexpr auto select_tags(
        UniformMesh<TagsHead, TagsQueue...> m,
        SliceSpecsHead&& h,
        SliceSpecsQueue&&... q) noexcept
{
    static_assert(sizeof...(TagsQueue) == sizeof...(SliceSpecsQueue));
    if constexpr (sizeof...(TagsQueue) > 0) {
        return append_if_all(
                select_tags(UniformMesh<TagsQueue...>(m), q...),
                UniformMesh<TagsHead>(m),
                m,
                h);
    } else {
        return append_if_all(UniformMesh<>(RCoord<>(), RLength<>()), m, m, h);
    }
}

} // namespace detail

template <class... Tags>
std::ostream& operator<<(std::ostream& out, UniformMesh<Tags...> const& dom)
{
    return out << "UniformMesh( origin=" << dom.origin() << ", unitvec=" << dom.step() << " )";
}

template <class... Tags, class... SliceSpecs>
inline constexpr auto submesh(const UniformMesh<Tags...>& mesh, SliceSpecs... slices) noexcept
{
    using ReturnType = decltype(detail::select_tags(mesh, std::forward<SliceSpecs>(slices)...));
    return ReturnType(mesh);
}

template <class... Tags>
using Mesh = UniformMesh<Tags...>;

using MeshX = Mesh<Dim::X>;

using MeshVx = Mesh<Dim::Vx>;

using MeshXVx = Mesh<Dim::X, Dim::Vx>;
