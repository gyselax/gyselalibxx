// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include <Kokkos_Core.hpp>

#include "ddc/coordinate.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_space.hpp"
#include "ddc/discrete_vector.hpp"
#include "ddc/real_type.hpp"

namespace ddc {

struct PeriodicSamplingBase
{
};

/** PeriodicSampling models a periodic discretization of the provided continuous dimension
 */
template <class CDim>
class PeriodicSampling : PeriodicSamplingBase
{
public:
    using continuous_dimension_type = CDim;

    using continuous_element_type = Coordinate<CDim>;


    using discrete_dimension_type = PeriodicSampling;

public:
    template <class DDim, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

    private:
        continuous_element_type m_origin {0};

        Real m_step {1};

        std::size_t m_n_period {2};

    public:
        using discrete_dimension_type = PeriodicSampling;

        using discrete_domain_type = DiscreteDomain<DDim>;

        using discrete_element_type = DiscreteElement<DDim>;

        using discrete_vector_type = DiscreteVector<DDim>;

        Impl() = default;

        Impl(Impl const&) = delete;

        template <class OriginMemorySpace>
        explicit Impl(Impl<DDim, OriginMemorySpace> const& impl)
            : m_origin(impl.m_origin)
            , m_step(impl.m_step)
            , m_n_period(impl.m_n_period)
        {
        }

        Impl(Impl&&) = default;

        /** @brief Construct a `Impl` from a point and a spacing step.
         *
         * @param origin the real coordinate of mesh coordinate 0
         * @param step   the real distance between two points of mesh distance 1
         * @param n_period   the number of steps in a period
         */
        Impl(continuous_element_type origin, Real step, std::size_t n_period)
            : m_origin(origin)
            , m_step(step)
            , m_n_period(n_period)
        {
            assert(step > 0);
            assert(n_period > 0);
        }

        ~Impl() = default;

        /// @brief Lower bound index of the mesh
        KOKKOS_FUNCTION continuous_element_type origin() const noexcept
        {
            return m_origin;
        }

        /// @brief Lower bound index of the mesh
        KOKKOS_FUNCTION discrete_element_type front() const noexcept
        {
            return discrete_element_type {0};
        }

        /// @brief Spacing step of the mesh
        KOKKOS_FUNCTION Real step() const
        {
            return m_step;
        }

        /// @brief Number of steps in a period
        KOKKOS_FUNCTION std::size_t n_period() const
        {
            return m_n_period;
        }

        /// @brief Convert a mesh index into a position in `CDim`
        KOKKOS_FUNCTION continuous_element_type
        coordinate(discrete_element_type const& icoord) const noexcept
        {
            return m_origin
                   + continuous_element_type(
                             static_cast<int>((icoord.uid() + m_n_period / 2 - 1) % m_n_period)
                             - static_cast<int>(m_n_period / 2 - 1))
                             * m_step;
        }
    };

    /** Construct a Impl<Kokkos::HostSpace> and associated discrete_domain_type from a segment
     *  \f$[a, b] \subset [a, +\infty[\f$ and a number of points `n`.
     *
     * @param a coordinate of the first point of the domain
     * @param b coordinate of the last point of the domain
     * @param n number of points to map on the segment \f$[a, b]\f$ including a & b
     * @param n_period   the number of steps in a period
     */
    template <class DDim>
    static std::tuple<typename DDim::template Impl<DDim, Kokkos::HostSpace>, DiscreteDomain<DDim>>
    init(Coordinate<CDim> a,
         Coordinate<CDim> b,
         DiscreteVector<DDim> n,
         DiscreteVector<DDim> n_period)
    {
        assert(a < b);
        assert(n > 1);
        assert(n_period > 1);
        typename DDim::template Impl<DDim, Kokkos::HostSpace>
                disc(a, Coordinate<CDim> {(b - a) / (n - 1)}, n_period);
        DiscreteDomain<DDim> domain {disc.front(), n};
        return std::make_tuple(std::move(disc), std::move(domain));
    }

    /** Construct a periodic `DiscreteDomain` from a segment \f$[a, b] \subset [a, +\infty[\f$ and a
     *  number of points `n`.
     *
     * @param a coordinate of the first point of the domain
     * @param b coordinate of the last point of the domain
     * @param n the number of points to map the segment \f$[a, b]\f$ including a & b
     * @param n_period   the number of steps in a period
     * @param n_ghosts_before number of additional "ghost" points before the segment
     * @param n_ghosts_after number of additional "ghost" points after the segment
     */
    template <class DDim>
    std::tuple<
            Impl<DDim, Kokkos::HostSpace>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>>
    init_ghosted(
            continuous_element_type a,
            continuous_element_type b,
            DiscreteVector<DDim> n,
            DiscreteVector<DDim> n_period,
            DiscreteVector<DDim> n_ghosts_before,
            DiscreteVector<DDim> n_ghosts_after)
    {
        using discrete_domain_type = DiscreteDomain<DDim>;
        assert(a < b);
        assert(n > 1);
        assert(n_period > 1);
        Real discretization_step {(b - a) / (n - 1)};
        Impl<DDim, Kokkos::HostSpace>
                disc(a - n_ghosts_before.value() * discretization_step,
                     discretization_step,
                     n_period);
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

    /** Construct a periodic `DiscreteDomain` from a segment \f$[a, b] \subset [a, +\infty[\f$ and a
     *  number of points `n`.
     *
     * @param a coordinate of the first point of the domain
     * @param b coordinate of the last point of the domain
     * @param n the number of points to map the segment \f$[a, b]\f$ including a & b
     * @param n_period   the number of steps in a period
     * @param n_ghosts number of additional "ghost" points before and after the segment
     */
    template <class DDim>
    std::tuple<
            Impl<DDim, Kokkos::HostSpace>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>,
            DiscreteDomain<DDim>>
    init_ghosted(
            continuous_element_type a,
            continuous_element_type b,
            DiscreteVector<DDim> n,
            DiscreteVector<DDim> n_period,
            DiscreteVector<DDim> n_ghosts)
    {
        return init_ghosted(a, b, n, n_period, n_ghosts, n_ghosts);
    }
};

template <class DDim>
struct is_periodic_sampling : public std::is_base_of<PeriodicSamplingBase, DDim>
{
};

template <class DDim>
constexpr bool is_periodic_sampling_v = is_periodic_sampling<DDim>::value;

template <
        class DDimImpl,
        std::enable_if_t<
                is_periodic_sampling_v<typename DDimImpl::discrete_dimension_type>,
                int> = 0>
std::ostream& operator<<(std::ostream& out, DDimImpl const& mesh)
{
    return out << "PeriodicSampling( origin=" << mesh.origin() << ", step=" << mesh.step() << " )";
}

/// @brief Lower bound index of the mesh
template <class DDim>
KOKKOS_FUNCTION std::enable_if_t<
        is_periodic_sampling_v<DDim>,
        Coordinate<typename DDim::continuous_dimension_type>>
origin() noexcept
{
    return discrete_space<DDim>().origin();
}

/// @brief Lower bound index of the mesh
template <class DDim>
KOKKOS_FUNCTION std::enable_if_t<is_periodic_sampling_v<DDim>, DiscreteElement<DDim>>
front() noexcept
{
    return discrete_space<DDim>().front();
}

/// @brief Spacing step of the mesh
template <class DDim>
KOKKOS_FUNCTION std::enable_if_t<is_periodic_sampling_v<DDim>, Real> step() noexcept
{
    return discrete_space<DDim>().step();
}

template <class DDim, std::enable_if_t<is_periodic_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> coordinate(
        DiscreteElement<DDim> const& c)
{
    return discrete_space<DDim>().coordinate(c);
}

template <class DDim, std::enable_if_t<is_periodic_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> distance_at_left(
        DiscreteElement<DDim>)
{
    return Coordinate<typename DDim::continuous_dimension_type>(step<DDim>());
}

template <class DDim, std::enable_if_t<is_periodic_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> distance_at_right(
        DiscreteElement<DDim>)
{
    return Coordinate<typename DDim::continuous_dimension_type>(step<DDim>());
}

template <class DDim, std::enable_if_t<is_periodic_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> rmin(
        DiscreteDomain<DDim> const& d)
{
    return coordinate(d.front());
}

template <class DDim, std::enable_if_t<is_periodic_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> rmax(
        DiscreteDomain<DDim> const& d)
{
    return coordinate(d.back());
}

template <class DDim, std::enable_if_t<is_periodic_sampling_v<DDim>, int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type> rlength(
        DiscreteDomain<DDim> const& d)
{
    return rmax(d) - rmin(d);
}

} // namespace ddc
