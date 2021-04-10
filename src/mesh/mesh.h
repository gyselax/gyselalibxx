#pragma once

#include <cassert>
#include <functional>
#include <initializer_list>
#include <tuple>

#include "space.h"
#include "view.h"

using MCoord = size_t;

template <int NDIMS>
using MCoordND = std::array<MCoord, NDIMS>;

using MCoord1D = MCoordND<1>;

using MCoord2D = MCoordND<2>;

using MCoord3D = MCoordND<3>;

class Mesh
{
    /// origin
    RCoord m_origin;

    /// step size
    double m_step;

public:
    constexpr Mesh(RCoord origin, double step) noexcept : m_origin(origin), m_step(step) { }

    constexpr double origin() const noexcept
    {
        return m_origin;
    }

    constexpr double step() const noexcept
    {
        return m_step;
    }

    RCoord1D operator()(const MCoord1D icoord) const
    {
        return {origin() + icoord[0] * m_step};
    }

    void operator()(View1D<RCoord1D> rcoords, const View1D<MCoord1D> icoords) const
    {
        assert(icoords.extents() == rcoords.extents());
        for (int ii = 0; ii < icoords.extent(0); ++ii) {
            rcoords[ii][0] = origin() + icoords[ii][0] * m_step;
        }
    }
};

template <int NDIMS>
using MeshND = std::array<Mesh, NDIMS>;

using Mesh1D = MeshND<1>;

using Mesh2D = MeshND<2>;

using Mesh3D = MeshND<3>;

class MDomain
{
public:
    using index_type = ptrdiff_t;

private:
    Mesh m_mesh;

    ptrdiff_t m_size;

public:
    /** Constructs a new ::MDomainND on the provided mesh.
     */
    constexpr MDomain(Mesh mesh, ptrdiff_t size) noexcept : m_mesh(mesh), m_size(size) { }

    constexpr MDomain(const MDomain&) noexcept = default;

    constexpr MDomain(MDomain&&) noexcept = default;

    constexpr MDomain& operator=(const MDomain&) noexcept = default;

    constexpr MDomain& operator=(MDomain&&) noexcept = default;

    constexpr const Mesh& mesh() const noexcept
    {
        return m_mesh;
    }

    constexpr void rcoords(View1D<RCoord1D> coords) const noexcept
    {
        for (size_t ii = 0; ii < coords.extent(0); ++ii) {
            coords[ii] = mesh()({ii});
        }
    }

    constexpr ptrdiff_t size() const noexcept
    {
        return m_size;
    }
};

template <size_t NDIM>
using MDomainND = std::array<MDomain, NDIM>;

using MDomain1D = MDomainND<1>;

using MDomain2D = MDomainND<2>;

using MDomain3D = MDomainND<3>;
