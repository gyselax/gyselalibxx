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
 * @param[in] domain
 * 	The domain on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the trapezoid method defined on the provided domain.
 */
template <class ExecSpace, class IDim>
device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<IDim>>> trapezoid_quadrature_coefficients_1d(
        ddc::DiscreteDomain<IDim> const& domain,
        device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<IDim>>> coefficients)
{
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
                coefficients(domain.front()) = 0.5 * dx_r; // distance_at_right(domain.front());
                coefficients(domain.back()) = 0.5 * dx_l; //distance_at_left(domain.back());
                if constexpr (IDim::continuous_dimension_type::PERIODIC) {
                    coefficients(domain.front()) += 0.5 * distance_at_left(domain.back());
                    coefficients(domain.back()) += 0.5 * distance_at_right(domain.front());
                }
            });
    return coefficients;
}

/**
 * @brief Get the trapezoid coefficients in ND.
 *
 * Calculate the quadrature coefficients for the trapezoid method defined on the provided domain.
 *
 * @param[in] domain
 * 	The domain on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the trapezoid method defined on the provided domain.
 */
template <class ExecSpace, class... ODims>
device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<ODims...>>> trapezoid_quadrature_coefficients(
        ddc::DiscreteDomain<ODims...> const& domain,
        device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<ODims...>>> coeffs)
{
    return quadrature_coeffs_nd<ExecSpace, ODims...>(
            domain,
            coeffs,
            (std::function<device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<ODims>>>(
                     ddc::DiscreteDomain<ODims>,
                     device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<ODims>>>)>(
                    trapezoid_quadrature_coefficients_1d<ExecSpace, ODims>))...);
}
