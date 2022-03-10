// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <vector>

#include "ddc/coordinate.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discretization.hpp"

/// `NonUniformDiscretization` models a non-uniform discretization of the `CDim` segment \f$[a, b]\f$.
template <class CDim>
class NonUniformDiscretization
{
public:
    using rcoord_type = Coordinate<CDim>;

    using mcoord_type = DiscreteCoordinate<NonUniformDiscretization>;

    using rdim_type = CDim;

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

    NonUniformDiscretization(NonUniformDiscretization const& x) = delete;

    NonUniformDiscretization(NonUniformDiscretization&& x) = delete;

    ~NonUniformDiscretization() = default;

    constexpr std::size_t size() const
    {
        return m_points.size();
    }

    /// @brief Convert a mesh index into a position in `CDim`
    constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        return m_points[icoord.value()];
    }
};

template <class>
struct is_non_uniform_discretization : public std::false_type
{
};

template <class CDim>
struct is_non_uniform_discretization<NonUniformDiscretization<CDim>> : public std::true_type
{
};

template <class DDim>
constexpr bool is_non_uniform_discretization_v = is_non_uniform_discretization<DDim>::value;

template <class CDim>
std::ostream& operator<<(std::ostream& out, NonUniformDiscretization<CDim> const& mesh)
{
    return out << "NonUniformDiscretization(" << mesh.size() << ")";
}

template <class CDim>
Coordinate<CDim> to_real(DiscreteCoordinate<NonUniformDiscretization<CDim>> const& c)
{
    return discretization<NonUniformDiscretization<CDim>>().to_real(c);
}

template <class CDim>
Coordinate<CDim> distance_at_left(DiscreteCoordinate<NonUniformDiscretization<CDim>> i)
{
    return to_real(i) - to_real(i - 1);
}

template <class CDim>
Coordinate<CDim> distance_at_right(DiscreteCoordinate<NonUniformDiscretization<CDim>> i)
{
    return to_real(i + 1) - to_real(i);
}

template <class CDim>
Coordinate<CDim> rmin(DiscreteDomain<NonUniformDiscretization<CDim>> const& d)
{
    return to_real(d.front());
}

template <class CDim>
Coordinate<CDim> rmax(DiscreteDomain<NonUniformDiscretization<CDim>> const& d)
{
    return to_real(d.back());
}

template <class CDim>
Coordinate<CDim> rlength(DiscreteDomain<NonUniformDiscretization<CDim>> const& d)
{
    return rmax(d) - rmin(d);
}
