#pragma once

#include <cassert>

#include "ddc/coordinate.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_dimension.hpp"

/** UniformDiscretization models a uniform discretization of the provided continuous dimension
 */
template <class CDim>
class UniformDiscretization : public DiscreteDimension
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
    rcoord_type m_origin = 0.;

    rcoord_type m_step = 1.;

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

    UniformDiscretization(UniformDiscretization const& x) = default;

    UniformDiscretization(UniformDiscretization&& x) = default;

    ~UniformDiscretization() = default;

    UniformDiscretization& operator=(UniformDiscretization const& x) = default;

    UniformDiscretization& operator=(UniformDiscretization&& x) = default;

    constexpr bool operator==(UniformDiscretization const& other) const
    {
        return m_origin == other.m_origin && m_step == other.m_step;
    }

    template <class ORDim>
    constexpr bool operator==(UniformDiscretization<ORDim> const& other) const
    {
        return false;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(UniformDiscretization const& other) const
    {
        return !(*this == other);
    }

    template <class ORDim>
    constexpr bool operator!=(UniformDiscretization<ORDim> const& other) const
    {
        return !(*this == other);
    }
#endif

    /// @brief Lower bound index of the mesh
    constexpr mcoord_type lbound() const noexcept
    {
        return mcoord_type(0);
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

    /// @brief Convert a mesh index into a position in `CDim`
    constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        assert(icoord >= lbound());
        return m_origin + rcoord_type(icoord.value()) * m_step;
    }

    /// @brief Position of the lower bound in `CDim`
    constexpr rcoord_type rmin() const noexcept
    {
        return m_origin;
    }
};

template <class CDim>
std::ostream& operator<<(std::ostream& out, UniformDiscretization<CDim> const& mesh)
{
    return out << "UniformDiscretization( origin=" << mesh.origin() << ", step=" << mesh.step()
               << " )";
}
