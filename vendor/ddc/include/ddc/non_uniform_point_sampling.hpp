// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

#include "ddc/coordinate.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_space.hpp"
#include "ddc/discrete_vector.hpp"

namespace ddc {

namespace detail {

struct NonUniformPointSamplingBase
{
};

} // namespace detail

/// `NonUniformPointSampling` models a non-uniform discretization of the `CDim` segment \f$[a, b]\f$.
template <class CDim>
class NonUniformPointSampling : detail::NonUniformPointSamplingBase
{
public:
    using continuous_dimension_type = CDim;

    using continuous_element_type = Coordinate<CDim>;


    using discrete_dimension_type = NonUniformPointSampling;

public:
    template <class DDim, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

        Kokkos::View<continuous_element_type*, MemorySpace> m_points;

    public:
        using discrete_dimension_type = NonUniformPointSampling;

        using discrete_domain_type = DiscreteDomain<DDim>;

        using discrete_element_type = DiscreteElement<DDim>;

        using discrete_vector_type = DiscreteVector<DDim>;

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
        explicit Impl(InputRange const& points)
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
        Impl(InputIt points_begin, InputIt points_end)
        {
            std::vector<continuous_element_type> host_points(points_begin, points_end);
            Kokkos::View<continuous_element_type*, Kokkos::HostSpace>
                    host(host_points.data(), host_points.size());
            Kokkos::resize(m_points, host.extent(0));
            Kokkos::deep_copy(m_points, host);
        }

        template <class OriginMemorySpace>
        explicit Impl(Impl<DDim, OriginMemorySpace> const& impl)
            : m_points(Kokkos::create_mirror_view_and_copy(MemorySpace(), impl.m_points))
        {
        }

        Impl(Impl const& x) = delete;

        Impl(Impl&& x) = default;

        ~Impl() = default;

        KOKKOS_FUNCTION std::size_t size() const
        {
            return m_points.size();
        }

        /// @brief Lower bound index of the mesh
        KOKKOS_FUNCTION discrete_element_type front() const noexcept
        {
            return discrete_element_type {0};
        }

        /// @brief Convert a mesh index into a position in `CDim`
        KOKKOS_FUNCTION continuous_element_type
        coordinate(discrete_element_type const& icoord) const noexcept
        {
            return m_points(icoord.uid());
        }
    };

    /** Construct an Impl<Kokkos::HostSpace> and associated discrete_domain_type from an iterator
     * containing the points coordinates along the `DDim` dimension.
     *
     * @param non_uniform_points a vector containing the coordinates of the points of the domain.
     */
    template <class DDim, class InputRange>
    static std::tuple<typename DDim::template Impl<DDim, Kokkos::HostSpace>, DiscreteDomain<DDim>>
    init(InputRange const non_uniform_points)
    {
        auto a = non_uniform_points.begin();
        auto b = non_uniform_points.end();
        auto n = std::distance(non_uniform_points.begin(), non_uniform_points.end());
        assert(a < b);
        assert(n > 0);
        typename DDim::template Impl<DDim, Kokkos::HostSpace> disc(non_uniform_points);
        DiscreteDomain<DDim> domain {disc.front(), DiscreteVector<DDim> {n}};
        return std::make_tuple(std::move(disc), std::move(domain));
    }

    /** Construct 4 non-uniform `DiscreteDomain` and an Impl<Kokkos::HostSpace> from 3 iterators containing the points coordinates along the `DDim` dimension.
     *
     * @param domain_r an iterator containing the coordinates of the points of the main domain along the DDim position
     * @param pre_ghost_r an iterator containing the positions of the ghost points before the main domain the DDim position
     * @param post_ghost_r an iterator containing the positions of the ghost points after the main domain the DDim position
     */
    template <class DDim, class InputRange>
    static std::tuple<
            typename DDim::template Impl<DDim, Kokkos::HostSpace>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>>
    init_ghosted(
            InputRange const& domain_r,
            InputRange const& pre_ghost_r,
            InputRange const& post_ghost_r)
    {
        using discrete_domain_type = DiscreteDomain<DDim>;
        auto n = DiscreteVector<DDim> {std::distance(domain_r.begin(), domain_r.end())};

        assert(domain_r.begin() < domain_r.end());
        assert(n > 1);

        auto n_ghosts_before
                = DiscreteVector<DDim> {std::distance(pre_ghost_r.begin(), pre_ghost_r.end())};
        auto n_ghosts_after
                = DiscreteVector<DDim> {std::distance(post_ghost_r.begin(), post_ghost_r.end())};

        std::vector<typename InputRange::value_type> full_domain;

        std::copy(pre_ghost_r.begin(), pre_ghost_r.end(), std::back_inserter(full_domain));
        std::copy(domain_r.begin(), domain_r.end(), std::back_inserter(full_domain));
        std::copy(post_ghost_r.begin(), post_ghost_r.end(), std::back_inserter(full_domain));

        typename DDim::template Impl<DDim, Kokkos::HostSpace> disc(full_domain);

        discrete_domain_type ghosted_domain
                = discrete_domain_type(disc.front(), n + n_ghosts_before + n_ghosts_after);
        discrete_domain_type pre_ghost
                = discrete_domain_type(ghosted_domain.front(), n_ghosts_before);
        discrete_domain_type main_domain
                = discrete_domain_type(ghosted_domain.front() + n_ghosts_before, n);
        discrete_domain_type post_ghost
                = discrete_domain_type(main_domain.back() + 1, n_ghosts_after);
        return std::make_tuple(
                std::move(disc),
                std::move(main_domain),
                std::move(ghosted_domain),
                std::move(pre_ghost),
                std::move(post_ghost));
    }
};

template <class DDim>
struct is_non_uniform_point_sampling
    : public std::is_base_of<detail::NonUniformPointSamplingBase, DDim>
{
};

template <class DDim>
constexpr bool is_non_uniform_point_sampling_v = is_non_uniform_point_sampling<DDim>::value;

template <class DDim>
using is_non_uniform_sampling [[deprecated("Use is_non_uniform_point_sampling instead")]]
= is_non_uniform_point_sampling<DDim>;

template <class DDim>
[[deprecated(
        "Use is_non_uniform_point_sampling_v instead")]] constexpr bool is_non_uniform_sampling_v
        = is_non_uniform_point_sampling_v<DDim>;

template <
        class DDimImpl,
        std::enable_if_t<
                is_non_uniform_point_sampling_v<typename DDimImpl::discrete_dimension_type>,
                int> = 0>
std::ostream& operator<<(std::ostream& out, DDimImpl const& mesh)
{
    return out << "NonUniformPointSampling(" << mesh.size() << ")";
}

template <class DDim, std::enable_if_t<is_non_uniform_point_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> coordinate(
        DiscreteElement<DDim> const& c)
{
    return discrete_space<DDim>().coordinate(c);
}

template <class DDim, std::enable_if_t<is_non_uniform_point_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> distance_at_left(
        DiscreteElement<DDim> i)
{
    return coordinate(i) - coordinate(i - 1);
}

template <class DDim, std::enable_if_t<is_non_uniform_point_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> distance_at_right(
        DiscreteElement<DDim> i)
{
    return coordinate(i + 1) - coordinate(i);
}

template <class DDim, std::enable_if_t<is_non_uniform_point_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> rmin(
        DiscreteDomain<DDim> const& d)
{
    return coordinate(d.front());
}

template <class DDim, std::enable_if_t<is_non_uniform_point_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> rmax(
        DiscreteDomain<DDim> const& d)
{
    return coordinate(d.back());
}

template <class DDim, std::enable_if_t<is_non_uniform_point_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> rlength(
        DiscreteDomain<DDim> const& d)
{
    return rmax(d) - rmin(d);
}

} // namespace ddc
