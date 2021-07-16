#pragma once

#include <cassert>

#include "mcoord.h"
#include "mesh.h"
#include "rcoord.h"
#include "taggedarray.h"

/// @class UniformMesh
/// @brief models a uniform discretization of the `RDim` half-line \f$[O, +\infty[\f$.
template <class RDim>
class UniformMesh : public Mesh
{
public:
    using rcoord_type = RCoord<RDim>;

    using rlength_type = RLength<RDim>;

    using mcoord_type = MCoord<UniformMesh>;

    using rdim_type = RDim;

public:
    static constexpr std::size_t rank()
    {
        return 1;
    }

private:
    double m_origin = 0.;

    double m_step = 1.;

public:
    UniformMesh() = default;

    /// @brief Construct a `UniformMesh` from a point and a spacing step.
    constexpr UniformMesh(rcoord_type origin, rcoord_type step) : m_origin(origin), m_step(step) {}

    /// @brief Construct a `UniformMesh` from a segment \f$[a, b] \subset [a, +\infty[\f$ and a number of points `n`.
    constexpr UniformMesh(rcoord_type a, rcoord_type b, std::size_t n)
        : m_origin(a)
        , m_step((b - a) / (n - 1))
    {
        assert(a < b);
        assert(n > 1);
    }

    UniformMesh(UniformMesh const& x) = default;

    UniformMesh(UniformMesh&& x) = default;

    ~UniformMesh() = default;

    UniformMesh& operator=(UniformMesh const& x) = default;

    UniformMesh& operator=(UniformMesh&& x) = default;

    constexpr bool operator==(UniformMesh const& other) const
    {
        return m_origin == other.m_origin && m_step == other.m_step;
    }

    template <class ORDim>
    constexpr bool operator==(UniformMesh<ORDim> const& other) const
    {
        return false;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(UniformMesh const& other) const
    {
        return !(*this == other);
    }

    template <class ORDim>
    constexpr bool operator!=(UniformMesh<ORDim> const& other) const
    {
        return !(*this == other);
    }
#endif

    /// @brief Lower bound index of the mesh
    constexpr mcoord_type lbound() const noexcept
    {
        return 0;
    }

    /// @brief Lower bound index of the mesh
    constexpr rcoord_type origin() const noexcept
    {
        return m_origin;
    }

    /// @brief Spacing step of the mesh
    constexpr rcoord_type step() const
    {
        return m_step;
    }

    /// @brief Convert a mesh index into a position in `RDim`
    constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        assert(icoord >= lbound());
        return m_origin + icoord * m_step;
    }

    /// @brief Position of the lower bound in `RDim`
    constexpr rcoord_type rmin() const noexcept
    {
        return m_origin;
    }
};

template <class RDim>
std::ostream& operator<<(std::ostream& out, UniformMesh<RDim> const& mesh)
{
    return out << "UniformMesh( origin=" << mesh.origin() << ", unitvec=" << mesh.step() << " )";
}
