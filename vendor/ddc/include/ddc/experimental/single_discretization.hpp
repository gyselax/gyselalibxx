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

    template <class MemorySpace>
    class Impl
    {
        template <class OMemorySpace>
        friend class Impl;

    private:
        /// origin
        rcoord_type m_point;

    public:
        inline constexpr Impl(rcoord_type origin) noexcept : m_point(std::move(origin)) {}

        Impl(Impl const& other) = delete;

        constexpr Impl(Impl&& other) = delete;

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
};

template <class Tag>
std::ostream& operator<<(std::ostream& out, SingleDiscretization<Tag> const& dom)
{
    return out << "SingleDiscretization( at=" << dom.origin() << " )";
}

} // namespace experimental

template <class>
struct is_single_discretization : public std::false_type
{
};

template <class CDim>
struct is_single_discretization<experimental::SingleDiscretization<CDim>> : public std::true_type
{
};

template <class DDim>
constexpr bool is_single_discretization_v = is_single_discretization<DDim>::value;


/// @brief coordinate of the mesh
template <class DDim>
std::enable_if_t<is_single_discretization_v<DDim>, typename DDim::rcoord_type> origin() noexcept
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return discretization_device<DDim>().origin();
#else
    return discretization_host<DDim>().origin();
#endif
}

template <class CDim>
Coordinate<CDim> to_real(DiscreteCoordinate<experimental::SingleDiscretization<CDim>> const& c)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return discretization_device<experimental::SingleDiscretization<CDim>>().to_real(c);
#else
    return discretization_host<experimental::SingleDiscretization<CDim>>().to_real(c);
#endif
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
