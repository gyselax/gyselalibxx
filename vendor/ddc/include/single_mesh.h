#pragma once

#include "mcoord.h"
#include "mesh.h"
#include "rcoord.h"
#include "taggedarray.h"

/// Experimental concept representing a mesh with a single point.  The
/// main purpose is when taking a subdomain in a multidimensional
/// mesh in order to loose track of the real position.
template <class RDim>
class SingleMesh : public Mesh
{
public:
    using rcoord_type = RCoord<RDim>;

    using mcoord_type = MCoord<SingleMesh>;

    using rdim_type = RDim;

private:
    /// origin
    rcoord_type m_point;

public:
    inline constexpr SingleMesh(SingleMesh const& other) noexcept : m_point(other.m_point) {}

    inline constexpr SingleMesh(SingleMesh&& other) noexcept : m_point(std::move(other.m_point)) {}

    inline constexpr SingleMesh(rcoord_type origin) noexcept : m_point(std::move(origin)) {}

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

    inline constexpr rcoord_type origin() const noexcept
    {
        return m_point;
    }

    static constexpr std::size_t size() noexcept
    {
        return 1;
    }

    inline constexpr rcoord_type to_real(const mcoord_type icoord) const noexcept
    {
        assert(icoord == mcoord_type(0));
        return m_point;
    }
};

template <class Tag>
std::ostream& operator<<(std::ostream& out, SingleMesh<Tag> const& dom)
{
    return out << "SinglePointMesh( origin=" << dom.origin() << " )";
}
