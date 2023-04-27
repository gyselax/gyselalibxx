// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

#include "ddc/coordinate.hpp"
#include "ddc/detail/macros.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_space.hpp"
#include "ddc/discrete_vector.hpp"

namespace ddc {

/// `NonUniformPointSampling` models a non-uniform discretization of the `CDim` segment \f$[a, b]\f$.
template <class CDim>
class NonUniformPointSampling
{
public:
    using continuous_dimension_type = CDim;

    using continuous_element_type = Coordinate<CDim>;


    using discrete_dimension_type = NonUniformPointSampling;

    using discrete_domain_type = DiscreteDomain<NonUniformPointSampling>;

    using discrete_element_type = DiscreteElement<NonUniformPointSampling>;

    using discrete_vector_type = DiscreteVector<NonUniformPointSampling>;

public:
    template <class MemorySpace>
    class Impl
    {
        template <class OMemorySpace>
        friend class Impl;

        Kokkos::View<continuous_element_type*, MemorySpace> m_points;

    public:
        using discrete_dimension_type = NonUniformPointSampling<CDim>;

        Impl() = default;

        /// @brief Construct a `NonUniformPointSampling` using a brace-list, i.e. `NonUniformPointSampling mesh({0., 1.})`
        explicit Impl(std::initializer_list<continuous_element_type> points)
        {
            std::vector<continuous_element_type> host_points(points.begin(), points.end());
            Kokkos::View<continuous_element_type*, Kokkos::HostSpace>
                    host(host_points.data(), host_points.size());
            Kokkos::resize(m_points, host.extent(0));
            Kokkos::deep_copy(m_points, host);
        }

        /// @brief Construct a `NonUniformPointSampling` using a C++20 "common range".
        template <class InputRange>
        explicit inline constexpr Impl(InputRange const& points)
        {
            if constexpr (Kokkos::is_view<InputRange>::value) {
                Kokkos::deep_copy(m_points, points);
            } else {
                std::vector<continuous_element_type> host_points(points.begin(), points.end());
                Kokkos::View<continuous_element_type*, Kokkos::HostSpace>
                        host(host_points.data(), host_points.size());
                Kokkos::resize(m_points, host.extent(0));
                Kokkos::deep_copy(m_points, host);
            }
        }

        /// @brief Construct a `NonUniformPointSampling` using a pair of iterators.
        template <class InputIt>
        inline constexpr Impl(InputIt points_begin, InputIt points_end)
        {
            std::vector<continuous_element_type> host_points(points_begin, points_end);
            Kokkos::View<continuous_element_type*, Kokkos::HostSpace>
                    host(host_points.data(), host_points.size());
            Kokkos::resize(m_points, host.extent(0));
            Kokkos::deep_copy(m_points, host);
        }

        template <class OriginMemorySpace>
        explicit Impl(Impl<OriginMemorySpace> const& impl)
            : m_points(Kokkos::create_mirror_view_and_copy(MemorySpace(), impl.m_points))
        {
        }

        Impl(Impl const& x) = delete;

        Impl(Impl&& x) = default;

        ~Impl() = default;

        constexpr std::size_t size() const
        {
            return m_points.size();
        }

        /// @brief Convert a mesh index into a position in `CDim`
        constexpr continuous_element_type coordinate(
                discrete_element_type const& icoord) const noexcept
        {
            return m_points(icoord.uid());
        }
    };
};

template <class>
struct is_non_uniform_sampling : public std::false_type
{
};

template <class CDim>
struct is_non_uniform_sampling<NonUniformPointSampling<CDim>> : public std::true_type
{
};

template <class DDim>
constexpr bool is_non_uniform_sampling_v = is_non_uniform_sampling<DDim>::value;

template <
        class DDimImpl,
        std::enable_if_t<
                is_non_uniform_sampling_v<typename DDimImpl::discrete_dimension_type>,
                int> = 0>
std::ostream& operator<<(std::ostream& out, DDimImpl const& mesh)
{
    return out << "NonUniformPointSampling(" << mesh.size() << ")";
}

template <class CDim>
DDC_INLINE_FUNCTION Coordinate<CDim> distance_at_left(
        DiscreteElement<NonUniformPointSampling<CDim>> i)
{
    return coordinate(i) - coordinate(i - 1);
}

template <class CDim>
DDC_INLINE_FUNCTION Coordinate<CDim> distance_at_right(
        DiscreteElement<NonUniformPointSampling<CDim>> i)
{
    return coordinate(i + 1) - coordinate(i);
}

template <class CDim>
DDC_INLINE_FUNCTION Coordinate<CDim> rmin(DiscreteDomain<NonUniformPointSampling<CDim>> const& d)
{
    return coordinate(d.front());
}

template <class CDim>
DDC_INLINE_FUNCTION Coordinate<CDim> rmax(DiscreteDomain<NonUniformPointSampling<CDim>> const& d)
{
    return coordinate(d.back());
}

template <class CDim>
DDC_INLINE_FUNCTION Coordinate<CDim> rlength(DiscreteDomain<NonUniformPointSampling<CDim>> const& d)
{
    return rmax(d) - rmin(d);
}

} // namespace ddc
