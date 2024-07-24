// SPDX-License-Identifier: MIT
/** @file trapezoid_quadrature.hpp
 * File providing quadrature coefficients via the trapezoidal method.
 */
#pragma once

#include <ddc/ddc.hpp>

#include "quadrature_coeffs_nd.hpp"

/**
 * @brief Get the trapezoid coefficients in 1D.
 *
 * Calculate the quadrature coefficients for the trapezoid method defined on the provided domain.
 *
 * @tparam ExecSpace Execution space, depends on Kokkos.
 *
 * @param[in] domain
 * 	The domain on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the trapezoid method defined on the provided domain.
 *         The allocation place (host or device ) will depend on the ExecSpace.
 */
template <class ExecSpace, class IDim>
ddc::Chunk<
        double,
        ddc::DiscreteDomain<IDim>,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
trapezoid_quadrature_coefficients_1d(ddc::DiscreteDomain<IDim> const& domain)
{
    ddc::Chunk<
            double,
            ddc::DiscreteDomain<IDim>,
            ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients_alloc(domain);
    ddc::ChunkSpan const coefficients = coefficients_alloc.span_view();
    ddc::DiscreteDomain<IDim> middle_domain
            = domain.remove(ddc::DiscreteVector<IDim>(1), ddc::DiscreteVector<IDim>(1));

    ddc::parallel_for_each(
            ExecSpace(),
            middle_domain,
            KOKKOS_LAMBDA(ddc::DiscreteElement<IDim> const idx) {
                coefficients(idx) = 0.5 * (distance_at_left(idx) + distance_at_right(idx));
            });
    double const dx_l = distance_at_left(domain.back());
    double const dx_r = distance_at_right(domain.front());
    Kokkos::parallel_for(
            "bounds",
            Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int i) {
                coefficients(domain.front()) = 0.5 * dx_r;
                coefficients(domain.back()) = 0.5 * dx_l;
                if constexpr (IDim::continuous_dimension_type::PERIODIC) {
                    coefficients(domain.front()) += 0.5 * dx_l;
                    coefficients(domain.back()) += 0.5 * dx_r;
                }
            });

    return coefficients_alloc;
}

/**
 * @brief Get the trapezoid coefficients in ND.
 *
 * Calculate the quadrature coefficients for the trapezoid method defined on the provided domain.
 *
 * @tparam ExecSpace Execution space, depends on Kokkos.
 *
 * @param[in] domain
 * 	The domain on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the trapezoid method defined on the provided domain.
 *         The allocation place (host or device ) will depend on the ExecSpace.
 */
template <class ExecSpace, class... ODims>
ddc::Chunk<
        double,
        ddc::DiscreteDomain<ODims...>,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
trapezoid_quadrature_coefficients(ddc::DiscreteDomain<ODims...> const& domain)
{
    return quadrature_coeffs_nd<ExecSpace, ODims...>(
            domain,
            (std::function<ddc::Chunk<
                     double,
                     ddc::DiscreteDomain<ODims>,
                     ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>(
                     ddc::DiscreteDomain<ODims>)>(
                    trapezoid_quadrature_coefficients_1d<ExecSpace, ODims>))...);
}