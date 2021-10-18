#pragma once

#include "ddc/coordinate.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_dimension.hpp"

namespace experimental {

/** Experimental concept representing a discretization with a single point.
 *
 * The main purpose is when taking a subdomain in a multidimensional discrete space in order to
 * loose track of the real position.
 */
template <class RDim>
class SingleDiscretization : public DiscreteDimension
{
public:
    using rcoord_type = Coordinate<RDim>;

    using mcoord_type = DiscreteCoordinate<SingleDiscretization>;

    using rdim_type = RDim;

private:
    /// origin
    rcoord_type m_point;

public:
    inline constexpr SingleDiscretization(SingleDiscretization const& other) noexcept
        : m_point(other.m_point)
    {
    }

    inline constexpr SingleDiscretization(SingleDiscretization&& other) noexcept
        : m_point(std::move(other.m_point))
    {
    }

    inline constexpr SingleDiscretization(rcoord_type origin) noexcept : m_point(std::move(origin))
    {
    }

    friend constexpr bool operator==(SingleDiscretization const& xx, SingleDiscretization const& yy)
    {
        return (&xx == &yy) || (xx.m_point == yy.m_point);
    }

    friend constexpr bool operator!=(SingleDiscretization const& xx, SingleDiscretization const& yy)
    {
        return (&xx != &yy) && (xx.m_point != yy.m_point);
    }

    template <class OTag>
    friend constexpr bool operator==(
            SingleDiscretization const& xx,
            SingleDiscretization<OTag> const& yy)
    {
        return false;
    }

    template <class OTag>
    friend constexpr bool operator!=(
            SingleDiscretization const& xx,
            SingleDiscretization<OTag> const& yy)
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
std::ostream& operator<<(std::ostream& out, SingleDiscretization<Tag> const& dom)
{
    return out << "SingleDiscretization( at=" << dom.origin() << " )";
}

} // namespace experimental
