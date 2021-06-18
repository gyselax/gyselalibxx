#pragma once

#include "mcoord.h"
#include "rcoord.h"
#include "taggedarray.h"

/// Experimental concept representing a mesh with a single point.  The
/// main purpose is when taking a subdomain in a multidimensional
/// mesh in order to loose track of the real position.
template <class Tag>
class SingleMesh
{
public:
    using RCoord_ = RCoord<Tag>;

    using MCoord_ = MCoord<Tag>;

    using Tag_ = Tag;

private:
    /// origin
    RCoord_ m_point;

public:
    inline constexpr SingleMesh(SingleMesh const& other) noexcept : m_point(other.m_point) {}

    inline constexpr SingleMesh(SingleMesh&& other) noexcept : m_point(std::move(other.m_point)) {}

    inline constexpr SingleMesh(RCoord_ origin) noexcept : m_point(std::move(origin)) {}

    friend constexpr bool operator==(SingleMesh const& xx, SingleMesh const& yy)
    {
        return (&xx == &yy) || (xx.m_point == yy.m_point);
    }

    friend constexpr bool operator!=(SingleMesh const& xx, SingleMesh const& yy)
    {
        return (&xx != &yy) && (xx.m_point != yy.m_point);
    }

    template <class OTag>
    friend constexpr bool operator==(SingleMesh const& xx, SingleMesh<OTag> const& yy)
    {
        return false;
    }

    template <class OTag>
    friend constexpr bool operator!=(SingleMesh const& xx, SingleMesh<OTag> const& yy)
    {
        return false;
    }

    static inline constexpr size_t rank() noexcept
    {
        return 0;
    }

    inline constexpr RCoord_ origin() const noexcept
    {
        return m_point;
    }

    static constexpr std::size_t size() noexcept
    {
        return 1;
    }

    inline constexpr RCoord_ to_real(const MCoord_ icoord) const noexcept
    {
        assert(icoord == MCoord_(0));
        return m_point;
    }
};

template <class Tag>
std::ostream& operator<<(std::ostream& out, SingleMesh<Tag> const& dom)
{
    return out << "SinglePointMesh( origin=" << dom.origin() << " )";
}
