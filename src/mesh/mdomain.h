#pragma once

#include "mcoord.h"
#include "rcoord.h"

template <class... Tags>
class RegularMesh
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
    friend class RegularMesh;

public:
    inline constexpr RegularMesh(RCoord_ origin, RLength_ step) noexcept
        : m_origin(origin)
        , m_step(step)
    {
    }

    template <class... OTags>
    inline constexpr RegularMesh(const RegularMesh<OTags...>& other) noexcept
        : m_origin(other.m_origin)
        , m_step(other.m_step)
    {
    }

    template <class... OTags>
    inline constexpr RegularMesh(RegularMesh<OTags...>&& other) noexcept
        : m_origin(std::move(other.m_origin))
        , m_step(std::move(other.m_step))
    {
    }

    inline constexpr RCoord_ origin() const noexcept
    {
        return m_origin;
    }

    inline constexpr RLength_ step() const noexcept
    {
        return m_step;
    }

    inline constexpr RCoord_ rcoord(const MCoord_ icoord) const noexcept
    {
        return {origin() + icoord * m_step};
    }

    inline constexpr RCoord_ operator()(const MCoord_ icoord) const noexcept
    {
        return rcoord(icoord);
    }

    static inline constexpr size_t rank() noexcept
    {
        return sizeof...(Tags);
    }
};


namespace detail {

template <class... SelectedTags, class TagsHead, class... TagsQueue, class SliceSpec>
inline constexpr auto append_if_all(RegularMesh<SelectedTags...>, RegularMesh<TagsHead>, SliceSpec ) noexcept
{
	static_assert(std::is_integral_v<SliceSpec> || std::is_same_v<std::experimental::all_type, SliceSpec> );
    if constexpr (std::is_integral_v<SliceSpec>) {
        return RegularMesh<SelectedTags...>();
    } else {
        return RegularMesh<TagsHead, SelectedTags...>();
    }
}


template <class TagsHead, class... TagsQueue, class SliceSpecsHead, class... SliceSpecsQueue>
inline constexpr auto select_tags(RegularMesh<TagsHead, TagsQueue...>, SliceSpecsHead h, SliceSpecsQueue... q) noexcept
{
	return append_if_all(select_tags(RegularMesh<TagsQueue...>(), q...), RegularMesh<TagsHead>(), h );
}

template <class Tag, class SliceSpec>
inline constexpr auto select_tags(RegularMesh<Tag> m, SliceSpec s) noexcept
{
	return append_if_all(RegularMesh<>(RCoord<>(), RLength<>()), m, s );
}

} // namespace detail

template <class... Tags, class... SliceSpecs>
inline constexpr auto submesh(const RegularMesh<Tags...>& mesh, SliceSpecs... slices) noexcept
{
	using ReturnType = decltype(detail::select_tags(mesh, std::forward<SliceSpecs>(slices)...));
    return ReturnType(mesh);
}


template <class... Tags>
class RegularMDomain : public RegularMesh<Tags...>
{
public:
    using RegularMesh_ = RegularMesh<Tags...>;

    using RCoord_ = RCoord<Tags...>;

    using MCoord_ = MCoord<Tags...>;

    using Mesh = RegularMesh_;

private:
    /// step size
    MCoord_ m_begin;

    /// step size
    MCoord_ m_end;

public:
    inline constexpr RegularMDomain(RegularMesh_ mesh, RCoord_ begin, RCoord_ end) noexcept
        : RegularMesh_(std::move(mesh))
        , m_begin(std::move(begin))
        , m_end(std::move(end))
    {
    }

    inline constexpr RegularMDomain(RegularMesh_ mesh, RCoord_ end) noexcept
        : RegularMesh_(std::move(mesh))
        , m_begin(0)
        , m_end(std::move(end))
    {
    }

    inline constexpr RegularMDomain(
            RCoord_ origin,
            double step,
            RCoord_ begin,
            RCoord_ end) noexcept
        : RegularMesh_(std::move(origin), step)
        , m_begin(std::move(begin))
        , m_end(std::move(end))
    {
    }

    inline constexpr RegularMDomain(RCoord_ origin, double step, RCoord_ end) noexcept
        : RegularMesh_(std::move(origin), step)
        , m_begin(0)
        , m_end(std::move(end))
    {
    }

    template <class... OTags>
    inline constexpr RegularMDomain(const RegularMDomain<OTags...>& other) noexcept
        : RegularMesh_(other)
        , m_begin(other.m_begin)
        , m_end(other.m_end)
    {
    }

    template <class... OTags>
    inline constexpr RegularMDomain(RegularMDomain<OTags...>&& other) noexcept
        : RegularMesh_(std::move(other))
        , m_begin(std::move(other.m_begin))
        , m_end(std::move(other.m_end))
    {
    }

    inline constexpr MCoord_& begin() noexcept
    {
        return m_begin;
    }

    inline constexpr const MCoord_& begin() const noexcept
    {
        return m_begin;
    }

    inline constexpr const MCoord_& cbegin() const noexcept
    {
        return m_begin;
    }

    inline constexpr MCoord_& end() noexcept
    {
        return m_end;
    }

    inline constexpr const MCoord_& end() const noexcept
    {
        return m_end;
    }

    inline constexpr const MCoord_& cend() const noexcept
    {
        return m_end;
    }

    template <class QueryTag>
    inline constexpr size_t extent() const noexcept
    {
        return get<QueryTag>(m_end) - get<QueryTag>(m_begin);
    }

    inline constexpr size_t size() const noexcept
    {
        return ((extent<Tags>()) * ...);
    }
};


/* For now MDomain is just an alias to RegularMDomain, in the long run, we should use a tuple-based
 * solutions to have different types in each dimension
 */
template <class... Tags>
using MDomain = RegularMDomain<Tags...>;

using MDomainX = MDomain<Dim::X>;

using MDomainVx = MDomain<Dim::Vx>;

using MDomainXVx = MDomain<Dim::X, Dim::Vx>;
