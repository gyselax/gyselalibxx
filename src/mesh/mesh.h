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
    inline constexpr Mesh(RCoord origin, double step) noexcept : m_origin(origin), m_step(step) { }

    inline constexpr double origin() const noexcept
    {
        return m_origin;
    }

    inline constexpr double step() const noexcept
    {
        return m_step;
    }

    inline constexpr RCoord1D operator()(const MCoord1D icoord) const noexcept
    {
        return {origin() + icoord[0] * m_step};
    }

    inline constexpr void operator()(View1D<RCoord1D> rcoords, const View1D<MCoord1D> icoords)
            const noexcept
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

namespace detail {

template <class... SliceSpecs>
inline constexpr size_t mesh_size() noexcept
{
    return ((std::is_same_v<SliceSpecs, std::experimental::all_type> ? 1 : 0) + ... + 0);
}

template <size_t... As, size_t... Bs>
inline constexpr std::index_sequence<As..., Bs...> cat(
        std::index_sequence<As...>,
        std::index_sequence<Bs...>) noexcept
{
    return {};
}

template <size_t IDX, class SliceSpec>
inline constexpr auto if_nondouble(std::index_sequence<IDX>, SliceSpec) noexcept
{
    if constexpr (std::is_integral_v<SliceSpec>) {
        return std::index_sequence<> {};
    } else {
        return std::index_sequence<IDX> {};
    }
}

inline constexpr std::index_sequence<> nondouble_idxs(std::index_sequence<>) noexcept
{
    return {};
}

template <size_t IDX0, size_t... IDXS, class SliceSpec0, class... SliceSpecs>
inline constexpr auto nondouble_idxs(
        std::index_sequence<IDX0, IDXS...>,
        SliceSpec0,
        SliceSpecs... slices) noexcept
{
    return cat(
            if_nondouble(std::index_sequence<IDX0> {}, SliceSpec0 {}),
            nondouble_idxs(std::index_sequence<IDXS...> {}, slices...));
}

template <size_t NDIMS, size_t... IDXS>
inline constexpr auto submesh_impl(const MeshND<NDIMS>& mesh, std::index_sequence<IDXS...>) noexcept
{
    return MeshND<sizeof...(IDXS)> {mesh[IDXS]...};
}

} // namespace detail

template <size_t NDIMS, class... SliceSpecs>
inline constexpr auto submesh(const MeshND<NDIMS>& mesh, SliceSpecs... slices) noexcept
{
    return detail::submesh_impl(
            mesh,
            detail::nondouble_idxs(
                    std::make_index_sequence<sizeof...(slices)>(),
                    std::forward<SliceSpecs>(slices)...));
}



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

    inline constexpr MDomain(Mesh mesh, ptrdiff_t size) noexcept : m_mesh(mesh), m_size(size) { }


    inline constexpr MDomain(const MDomain&) noexcept = default;


    inline constexpr MDomain(MDomain&&) noexcept = default;


    inline constexpr MDomain& operator=(const MDomain&) noexcept = default;


    inline constexpr MDomain& operator=(MDomain&&) noexcept = default;


    inline constexpr const Mesh& mesh() const noexcept
    {
        return m_mesh;
    }


    inline constexpr void rcoords(View1D<RCoord1D> coords) const noexcept
    {
        for (size_t ii = 0; ii < coords.extent(0); ++ii) {
            coords[ii] = mesh()({ii});
        }
    }


    inline constexpr ptrdiff_t size() const noexcept
    {
        return m_size;
    }
};

template <size_t NDIM>
using MDomainND = std::array<MDomain, NDIM>;

using MDomain1D = MDomainND<1>;

using MDomain2D = MDomainND<2>;

using MDomain3D = MDomainND<3>;
