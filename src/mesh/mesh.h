#pragma once

#include <cassert>
#include <functional>
#include <initializer_list>
#include <tuple>

#include "space.h"
#include "view.h"

template <int NDIMS>
using MCoordND = std::array<size_t, NDIMS>;

using MCoord1D = MCoordND<1>;

using MCoord2D = MCoordND<2>;

class Mesher
{
    /// origin coordinates
    double m_origin;

    /// step size
    double m_step;

public:
    Mesher(double origin, double step) : m_origin(origin), m_step(step) { }

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
        return {m_origin + icoord[0] * m_step};
    }

    void operator()(View1D<RCoord1D> rcoords, const View1D<MCoord1D> icoords) const
    {
        assert(icoords.extents() == rcoords.extents());
        for (int ii = 0; ii < icoords.extent(0); ++ii) {
            rcoords[ii][0] = m_origin + icoords[ii][0] * m_step;
        }
    }
};

template <int NDIMS>
using MeshND = std::array<Mesher, NDIMS>;

using Mesh1D = MeshND<1>;

using Mesh2D = MeshND<2>;

namespace detail {

template <class>
struct ExtentSize;

template <size_t RANK, size_t... ORANKS>
struct ExtentSize<std::index_sequence<RANK, ORANKS...>> {
    template <class Extents>
    static inline size_t eval(const Extents& e)
    {
        return e.extent(RANK) + ExtentSize<std::index_sequence<ORANKS...>>::eval(e);
    }
};

template <>
struct ExtentSize<std::index_sequence<>> {
    template <class Extents>
    static inline size_t eval(const Extents&)
    {
        return 0;
    }
};

template <class>
struct EndComputer;

template <size_t... RANKS>
struct EndComputer<std::index_sequence<RANKS...>> {
    template <class Extents>
    static inline MCoordND<sizeof...(RANKS)> eval(
            const MCoordND<sizeof...(RANKS)>& begin,
            const Extents& extents)
    {
        return MCoordND<sizeof...(RANKS)> {begin[RANKS] + extents.extent(RANKS)...};
    }
};

} // namespace detail

template <int NDIMS>
class MDomainND
{
public:
    using index_type = ptrdiff_t;

private:
    MCoordND<NDIMS> m_begin;

    ExtentsND<NDIMS> m_extents;

    MeshND<NDIMS> m_mesh;

public:
    /** Constructs a new NDomainND on the provided mesh.
     * 
     * The domain includes begin but excludes end
     */
    constexpr MDomainND(
            MeshND<NDIMS> mesh,
            MCoordND<NDIMS> begin,
            ExtentsND<NDIMS> extents) noexcept
        : m_begin(begin)
        , m_extents(extents)
        , m_mesh(mesh)
    {
    }

    constexpr MDomainND(const MDomainND&) noexcept = default;

    constexpr MDomainND(MDomainND&&) noexcept = default;

    constexpr MDomainND& operator=(const MDomainND&) noexcept = default;

    constexpr MDomainND& operator=(MDomainND&&) noexcept = default;

    constexpr const MeshND<NDIMS>& mesh() const noexcept
    {
        return m_mesh;
    }

    constexpr const Mesher& mesher(size_t dim) const noexcept
    {
        return m_mesh[dim];
    }

    constexpr void rcoords(View1D<RCoord1D> coords, size_t dim) const noexcept
    {
        for (size_t ii = 0; ii < m_extents.extent(dim); ++ii) {
            coords[ii] = m_mesh[dim]({ii + m_begin[dim]});
        }
    }

    constexpr const MCoordND<NDIMS>& begin() const noexcept
    {
        return m_begin;
    }

    constexpr MCoordND<NDIMS> end() const noexcept
    {
        return detail::EndComputer<std::make_index_sequence<NDIMS>>::eval(m_begin, m_extents);
    }

    constexpr const ExtentsND<NDIMS>& extents() const noexcept
    {
        return m_extents;
    }

    template <size_t... DIMS>
    constexpr MDomainND<sizeof...(DIMS)> slice() const noexcept
    {
        return MDomainND<sizeof...(DIMS)>(
                {m_mesh[DIMS]...},
                {m_begin[DIMS]...},
                ExtentsND<sizeof...(DIMS)> {m_extents.extent(DIMS)...});
    }

    static constexpr size_t rank() noexcept
    {
        return NDIMS;
    }

    static constexpr size_t rank_dynamic() noexcept
    {
        return NDIMS;
    }

    static constexpr index_type static_extent(size_t) noexcept
    {
        return std::experimental::dynamic_extent;
    }

    constexpr index_type extent(size_t dim) const noexcept
    {
        return m_extents.extent(dim);
    }

    constexpr ptrdiff_t size() const noexcept
    {
        return detail::ExtentSize<std::make_index_sequence<NDIMS>>::eval(m_extents);
    }
};

using MDomain1D = MDomainND<1>;

using MDomain2D = MDomainND<2>;
