#pragma once

#include <cassert>
#include <vector>

#include "ddc/mcoord.h"
#include "ddc/mesh.h"
#include "ddc/rcoord.h"

/// `NonUniformMesh` models a non-uniform discretization of the `RDim` segment \f$[a, b]\f$.
template <class RDim>
class NonUniformMesh : public Mesh
{
public:
    using rcoord_type = RCoord<RDim>;

    using rlength_type = RLength<RDim>;

    using mcoord_type = MCoord<NonUniformMesh>;

    using rdim_type = RDim;

public:
    static constexpr std::size_t rank()
    {
        return 1;
    }

private:
    std::vector<double> m_points;

public:
    NonUniformMesh() = default;

    /// @brief Construct a `NonUniformMesh` using a brace-list, i.e. `NonUniformMesh mesh({0., 1.})`
    explicit NonUniformMesh(std::initializer_list<rcoord_type> points)
        : m_points(points.begin(), points.end())
    {
    }

    /// @brief Construct a `NonUniformMesh` using a C++20 "common range".
    template <class InputRange>
    inline constexpr NonUniformMesh(InputRange&& points) : m_points(points.begin(), points.end())
    {
    }

    /// @brief Construct a `NonUniformMesh` using a pair of iterators.
    template <class InputIt>
    inline constexpr NonUniformMesh(InputIt points_begin, InputIt points_end)
        : m_points(points_begin, points_end)
    {
    }

    NonUniformMesh(NonUniformMesh const& x) = default;

    NonUniformMesh(NonUniformMesh&& x) = default;

    ~NonUniformMesh() = default;

    NonUniformMesh& operator=(NonUniformMesh const& x) = default;

    NonUniformMesh& operator=(NonUniformMesh&& x) = default;

    constexpr bool operator==(NonUniformMesh const& other) const
    {
        return m_points == other.m_points;
    }

    template <class ORDim>
    constexpr bool operator==(NonUniformMesh<ORDim> const& other) const
    {
        return false;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(NonUniformMesh const& other) const
    {
        return !(*this == other);
    }

    template <class ORDim>
    constexpr bool operator!=(NonUniformMesh<ORDim> const& other) const
    {
        return !(*this == other);
    }
#endif

    constexpr std::size_t size() const
    {
        return m_points.size();
    }

    /// @brief Lower bound index of the mesh
    constexpr mcoord_type lbound() const noexcept
    {
        return 0;
    }

    /// @brief Upper bound index of the mesh
    constexpr mcoord_type ubound() const noexcept
    {
        return m_points.size() - 1;
    }

    /// @brief Convert a mesh index into a position in `RDim`
    constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        assert(icoord >= lbound());
        assert(icoord <= ubound());
        return m_points[static_cast<std::size_t>(icoord)];
    }

    /// @brief Position of the lower bound in `RDim`
    constexpr rcoord_type rmin() const noexcept
    {
        return m_points.front();
    }

    /// @brief Position of the upper bound in `RDim`
    constexpr rcoord_type rmax() const noexcept
    {
        return m_points.back();
    }

    /// @brief Length of `RDim`
    constexpr rcoord_type rlength() const
    {
        return m_points.back() - m_points.front();
    }
};

template <class RDim>
std::ostream& operator<<(std::ostream& out, NonUniformMesh<RDim> const& mesh)
{
    out << "NonUniformMesh( ";
    if (mesh.size() > 0) {
        out << mesh.rmin() << "..." << mesh.rmax();
    }
    out << " )";
    return out;
}
