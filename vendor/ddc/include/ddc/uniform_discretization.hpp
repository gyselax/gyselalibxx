// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>

#include "ddc/coordinate.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discretization.hpp"

/** UniformDiscretization models a uniform discretization of the provided continuous dimension
 */
template <class CDim>
class UniformDiscretization
{
public:
    using rcoord_type = Coordinate<CDim>;

    using mcoord_type = DiscreteCoordinate<UniformDiscretization>;

    using rdim_type = CDim;

public:
    static constexpr std::size_t rank()
    {
        return 1;
    }

private:
    rcoord_type m_origin {0.};

    rcoord_type m_step {1.};

public:
    UniformDiscretization() = default;

    /** @brief Construct a `UniformDiscretization` from a point and a spacing step.
     * 
     * @param origin the real coordinate of mesh coordinate 0
     * @param step   the real distance between two points of mesh distance 1
     */
    constexpr UniformDiscretization(rcoord_type origin, rcoord_type step)
        : m_origin(origin)
        , m_step(step)
    {
        assert(step > 0);
    }

    /** @brief Construct a `UniformDiscretization` from a segment \f$[a, b] \subset [a, +\infty[\f$ and a number of points `n`.
     * 
     * @param a the coordinate of a first real point (will have mesh coordinate 0)
     * @param b the coordinate of the second real point (will have mesh coordinate `n-1`)
     * @param n the number of points to map the segment \f$[a, b]\f$ including a & b
     */
    constexpr UniformDiscretization(rcoord_type a, rcoord_type b, std::size_t n)
        : m_origin(a)
        , m_step((b - a) / (n - 1))
    {
        assert(a < b);
        assert(n > 1);
    }

    UniformDiscretization(UniformDiscretization const& x) = delete;

    UniformDiscretization(UniformDiscretization&& x) = delete;

    ~UniformDiscretization() = default;

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

    /// @brief Convert a mesh index into a position in `CDim`
    constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        return m_origin + rcoord_type(icoord.value()) * m_step;
    }
};

template <class>
struct is_uniform_discretization : public std::false_type
{
};

template <class CDim>
struct is_uniform_discretization<UniformDiscretization<CDim>> : public std::true_type
{
};

template <class DDim>
constexpr bool is_uniform_discretization_v = is_uniform_discretization<DDim>::value;


template <class CDim>
std::ostream& operator<<(std::ostream& out, UniformDiscretization<CDim> const& mesh)
{
    return out << "UniformDiscretization( origin=" << mesh.origin() << ", step=" << mesh.step()
               << " )";
}

/// @brief Lower bound index of the mesh
template <class DDim>
std::enable_if_t<is_uniform_discretization_v<DDim>, typename DDim::rcoord_type> origin() noexcept
{
    return discretization<DDim>().origin();
}

/// @brief Spacing step of the mesh
template <class DDim>
std::enable_if_t<is_uniform_discretization_v<DDim>, typename DDim::rcoord_type> step() noexcept
{
    return discretization<DDim>().step();
}

template <class CDim>
Coordinate<CDim> to_real(DiscreteCoordinate<UniformDiscretization<CDim>> const& c)
{
    return discretization<UniformDiscretization<CDim>>().to_real(c);
}

template <class CDim>
Coordinate<CDim> distance_at_left(DiscreteCoordinate<UniformDiscretization<CDim>> i)
{
    return discretization<UniformDiscretization<CDim>>().step();
}

template <class CDim>
Coordinate<CDim> distance_at_right(DiscreteCoordinate<UniformDiscretization<CDim>> i)
{
    return discretization<UniformDiscretization<CDim>>().step();
}

template <class CDim>
Coordinate<CDim> rmin(DiscreteDomain<UniformDiscretization<CDim>> const& d)
{
    return to_real(d.front());
}

template <class CDim>
Coordinate<CDim> rmax(DiscreteDomain<UniformDiscretization<CDim>> const& d)
{
    return to_real(d.back());
}

template <class CDim>
Coordinate<CDim> rlength(DiscreteDomain<UniformDiscretization<CDim>> const& d)
{
    return rmax(d) - rmin(d);
}
