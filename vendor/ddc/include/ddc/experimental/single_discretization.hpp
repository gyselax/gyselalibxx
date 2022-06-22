// SPDX-License-Identifier: MIT

#pragma once

#include "ddc/coordinate.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_space.hpp"

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
    using continuous_dimension_type = CDim;

    using continuous_element_type = Coordinate<CDim>;


    using discrete_dimension_type = SingleDiscretization;

    using discrete_domain_type = DiscreteDomain<SingleDiscretization>;

    using discrete_element_type = DiscreteElement<SingleDiscretization>;

    using discrete_vector_type = DiscreteVector<SingleDiscretization>;

    template <class MemorySpace>
    class Impl
    {
        template <class OMemorySpace>
        friend class Impl;

    private:
        /// origin
        continuous_element_type m_point;

    public:
        using discrete_dimension_type = SingleDiscretization;

        inline constexpr Impl(continuous_element_type origin) noexcept : m_point(std::move(origin))
        {
        }

        Impl(Impl const& other) = delete;

        constexpr Impl(Impl&& other) = delete;

        inline constexpr continuous_element_type origin() const noexcept
        {
            return m_point;
        }

        inline constexpr continuous_element_type coordinate(
                [[maybe_unused]] discrete_element_type icoord) const noexcept
        {
            assert(icoord == discrete_element_type(0));
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
std::enable_if_t<is_single_discretization_v<DDim>, typename DDim::continuous_element_type>
origin() noexcept
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return discrete_space_device<DDim>().origin();
#else
    return discrete_space_host<DDim>().origin();
#endif
}

template <class CDim>
Coordinate<CDim> coordinate(DiscreteElement<experimental::SingleDiscretization<CDim>> const& c)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return discrete_space_device<experimental::SingleDiscretization<CDim>>().coordinate(c);
#else
    return discrete_space_host<experimental::SingleDiscretization<CDim>>().coordinate(c);
#endif
}

template <class CDim>
Coordinate<CDim> rmin(DiscreteDomain<experimental::SingleDiscretization<CDim>> const& d)
{
    return coordinate(d.front());
}

template <class CDim>
Coordinate<CDim> rmax(DiscreteDomain<experimental::SingleDiscretization<CDim>> const& d)
{
    return coordinate(d.back());
}

template <class CDim>
Coordinate<CDim> rlength(DiscreteDomain<experimental::SingleDiscretization<CDim>> const& d)
{
    return rmax(d) - rmin(d);
}
