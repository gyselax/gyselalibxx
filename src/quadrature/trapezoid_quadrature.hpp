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
template <class IDim>
ddc::Chunk<double, ddc::DiscreteDomain<IDim>> trapezoid_quadrature_coefficients_1d(
        ddc::DiscreteDomain<IDim> const& domain)
{
    ddc::Chunk<double, ddc::DiscreteDomain<IDim>> coefficients(domain);
    ddc::DiscreteDomain<IDim> middle_domain
            = domain.remove(ddc::DiscreteVector<IDim>(1), ddc::DiscreteVector<IDim>(1));

    coefficients(domain.front()) = 0.5 * distance_at_right(domain.front());
    ddc::for_each(middle_domain, [&](ddc::DiscreteElement<IDim> const idx) {
        coefficients(idx) = 0.5 * (distance_at_left(idx) + distance_at_right(idx));
    });
    coefficients(domain.back()) = 0.5 * distance_at_left(domain.back());

    if constexpr (IDim::continuous_dimension_type::PERIODIC) {
        coefficients(domain.front()) += 0.5 * distance_at_left(domain.back());
        coefficients(domain.back()) += 0.5 * distance_at_right(domain.front());
    }

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
template <class... ODims>
ddc::Chunk<double, ddc::DiscreteDomain<ODims...>> trapezoid_quadrature_coefficients(
        ddc::DiscreteDomain<ODims...> const& domain)
{
    return quadrature_coeffs_nd(
            domain,
            (std::function<ddc::Chunk<double, ddc::DiscreteDomain<ODims>>(
                     ddc::DiscreteDomain<ODims>)>(trapezoid_quadrature_coefficients_1d<ODims>))...);
}
