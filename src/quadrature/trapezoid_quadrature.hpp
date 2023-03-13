// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

template <class IDim>
ddc::Chunk<double, ddc::DiscreteDomain<IDim>> trapezoid_quadrature_coefficients(
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
