#pragma once

#include <cassert>
#include <vector>

#include "ddc/coordinate.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_dimension.hpp"

/// `NonUniformDiscretization` models a non-uniform discretization of the `RDim` segment \f$[a, b]\f$.
template <class RDim>
class NonUniformDiscretization : public DiscreteDimension
{
public:
    using rcoord_type = Coordinate<RDim>;

    using mcoord_type = DiscreteCoordinate<NonUniformDiscretization>;

    using rdim_type = RDim;

public:
    static constexpr std::size_t rank()
    {
        return 1;
    }

private:
    std::vector<rcoord_type> m_points;

public:
    NonUniformDiscretization() = default;

    /// @brief Construct a `NonUniformDiscretization` using a brace-list, i.e. `NonUniformDiscretization mesh({0., 1.})`
    explicit NonUniformDiscretization(std::initializer_list<rcoord_type> points)
        : m_points(points.begin(), points.end())
    {
    }

    /// @brief Construct a `NonUniformDiscretization` using a C++20 "common range".
    template <class InputRange>
    inline constexpr NonUniformDiscretization(InputRange&& points)
        : m_points(points.begin(), points.end())
    {
    }

    /// @brief Construct a `NonUniformDiscretization` using a pair of iterators.
    template <class InputIt>
    inline constexpr NonUniformDiscretization(InputIt points_begin, InputIt points_end)
        : m_points(points_begin, points_end)
    {
    }

    NonUniformDiscretization(NonUniformDiscretization const& x) = default;

    NonUniformDiscretization(NonUniformDiscretization&& x) = default;

    ~NonUniformDiscretization() = default;

    NonUniformDiscretization& operator=(NonUniformDiscretization const& x) = default;

    NonUniformDiscretization& operator=(NonUniformDiscretization&& x) = default;

    constexpr bool operator==(NonUniformDiscretization const& other) const
    {
        return m_points == other.m_points;
    }

    template <class ORDim>
    constexpr bool operator==(NonUniformDiscretization<ORDim> const& other) const
    {
        return false;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(NonUniformDiscretization const& other) const
    {
        return !(*this == other);
    }

    template <class ORDim>
    constexpr bool operator!=(NonUniformDiscretization<ORDim> const& other) const
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
        return mcoord_type(0);
    }

    /// @brief Upper bound index of the mesh
    constexpr mcoord_type ubound() const noexcept
    {
        return mcoord_type(m_points.size() - 1);
    }

    /// @brief Convert a mesh index into a position in `RDim`
    constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        assert(icoord >= lbound());
        assert(icoord <= ubound());
        return m_points[icoord.value()];
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
std::ostream& operator<<(std::ostream& out, NonUniformDiscretization<RDim> const& mesh)
{
    out << "NonUniformDiscretization( ";
    if (mesh.size() > 0) {
        out << mesh.rmin() << ", ..., " << mesh.rmax();
    }
    out << " )";
    return out;
}
