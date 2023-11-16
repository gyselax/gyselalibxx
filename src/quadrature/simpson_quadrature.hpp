// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>


/**
 * @brief Get the Simpson coefficients in 1D.
 *
 * Calculate the quadrature coefficients for the Simpson method defined on the provided domain.
 *
 * @param[in] domain
 * 	The domain on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the Simpson method defined on the provided domain.
 */
template <class IDim>
ddc::Chunk<double, ddc::DiscreteDomain<IDim>> simpson_quadrature_coefficients_1d(
        ddc::DiscreteDomain<IDim> const& domain)
{
    ddc::Chunk<double, ddc::DiscreteDomain<IDim>> coefficients(domain);
    ddc::DiscreteDomain<IDim> middle_domain
            = domain.remove(ddc::DiscreteVector<IDim>(1), ddc::DiscreteVector<IDim>(1));

    coefficients(domain.front()) = 1. / 3. * distance_at_right(domain.front());

    for (auto it = domain.begin() + 1; it < domain.end() - 1; it += 2) {
        ddc::DiscreteElement<IDim> idx = *it;
        coefficients(idx) = 2. / 3. * (distance_at_left(idx) + distance_at_right(idx));
        idx += ddc::DiscreteVector<IDim>(1);
        coefficients(idx) = 1. / 3. * (distance_at_left(idx) + distance_at_right(idx));
    }
    coefficients(domain.back()) = 1. / 3. * distance_at_left(domain.back());

    if constexpr (IDim::continuous_dimension_type::PERIODIC) {
        coefficients(domain.front()) += 2. / 3 * distance_at_left(domain.back());
        coefficients(domain.back()) += 2. / 3 * distance_at_right(domain.front());
    }


    return coefficients;
}
