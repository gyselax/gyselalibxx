// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief Get the Simpson coefficients in 1D.
 *
 * Calculate the quadrature coefficients for the Simpson method defined on the provided index range.
 *
 * @param[in] index range
 * 	The index range on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the Simpson method defined on the provided index range.
 */
template <class Grid>
host_t<FieldMem<double, IdxRange<Grid>>> simpson_quadrature_coefficients_1d(
        IdxRange<Grid> const& idx_range)
{
    host_t<FieldMem<double, IdxRange<Grid>>> coefficients(idx_range);

    coefficients(idx_range.front()) = 1. / 3. * distance_at_right(idx_range.front());

    for (auto it = idx_range.begin() + 1; it < idx_range.end() - 1; it += 2) {
        Idx<Grid> idx = *it;
        coefficients(idx) = 2. / 3. * (distance_at_left(idx) + distance_at_right(idx));
        idx += IdxStep<Grid>(1);
        coefficients(idx) = 1. / 3. * (distance_at_left(idx) + distance_at_right(idx));
    }
    coefficients(idx_range.back()) = 1. / 3. * distance_at_left(idx_range.back());

    if constexpr (Grid::continuous_dimension_type::PERIODIC) {
        coefficients(idx_range.front()) += 2. / 3 * distance_at_left(idx_range.back());
        coefficients(idx_range.back()) += 2. / 3 * distance_at_right(idx_range.front());
    }

    return coefficients;
}
