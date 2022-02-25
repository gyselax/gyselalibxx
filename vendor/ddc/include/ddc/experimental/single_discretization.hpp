// SPDX-License-Identifier: MIT

#pragma once

#include "ddc/coordinate.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discretization.hpp"

namespace experimental {

/** Experimental concept representing a discretization with a single point.
 *
 * The main purpose is when taking a subdomain in a multidimensional discrete space in order to
 * loose track of the real position.
 */
template <class CDim>
class SingleDiscretization
{
public:
    using rcoord_type = Coordinate<CDim>;

    using mcoord_type = DiscreteCoordinate<SingleDiscretization>;

    using rdim_type = CDim;

private:
    /// origin
    rcoord_type m_point;

public:
    inline constexpr SingleDiscretization(rcoord_type origin) noexcept : m_point(std::move(origin))
    {
    }

    SingleDiscretization(SingleDiscretization const& other) = delete;

    constexpr SingleDiscretization(SingleDiscretization&& other) = delete;

    inline constexpr rcoord_type origin() const noexcept
    {
        return m_point;
    }

    inline constexpr rcoord_type to_real([[maybe_unused]] mcoord_type icoord) const noexcept
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

template <class>
struct is_single_disretization : public std::false_type
{
};

template <class CDim>
struct is_single_disretization<experimental::SingleDiscretization<CDim>> : public std::true_type
{
};

template <class DDim>
constexpr bool is_single_disretization_v = is_single_disretization<DDim>::value;


/// @brief coordinate of the mesh
template <class DDim>
std::enable_if_t<is_single_disretization_v<DDim>, typename DDim::rcoord_type> origin() noexcept
{
    return discretization<DDim>().origin();
}

template <class CDim>
Coordinate<CDim> to_real(DiscreteCoordinate<experimental::SingleDiscretization<CDim>> const& c)
{
    return discretization<experimental::SingleDiscretization<CDim>>().to_real(c);
}

template <class CDim>
Coordinate<CDim> rmin(DiscreteDomain<experimental::SingleDiscretization<CDim>> const& d)
{
    return to_real(d.front());
}

template <class CDim>
Coordinate<CDim> rmax(DiscreteDomain<experimental::SingleDiscretization<CDim>> const& d)
{
    return to_real(d.back());
}

template <class CDim>
Coordinate<CDim> rlength(DiscreteDomain<experimental::SingleDiscretization<CDim>> const& d)
{
    return rmax(d) - rmin(d);
}
