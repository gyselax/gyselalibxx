// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

template <class IDim>
Chunk<double, DiscreteDomain<IDim>> trapezoid_quadrature_coefficients(
        DiscreteDomain<IDim> const& domain)
{
    Chunk<double, DiscreteDomain<IDim>> coefficients(domain);
    DiscreteDomain<IDim> middle_domain
            = domain.remove(DiscreteVector<IDim>(1), DiscreteVector<IDim>(1));

    coefficients(domain.front()) = 0.5 * distance_at_right(domain.front());
    for_each(middle_domain, [&](DiscreteElement<IDim> const idx) {
        coefficients(idx) = 0.5 * (distance_at_left(idx) + distance_at_right(idx));
    });
    coefficients(domain.back()) = 0.5 * distance_at_left(domain.back());

    if constexpr (IDim::continuous_dimension_type::PERIODIC) {
        coefficients(domain.front()) += 0.5 * distance_at_left(domain.back());
        coefficients(domain.back()) += 0.5 * distance_at_right(domain.front());
    }

    return coefficients;
}
