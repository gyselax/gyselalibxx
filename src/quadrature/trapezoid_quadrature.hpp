// SPDX-License-Identifier: MIT
/** @file trapezoid_quadrature.hpp
 * File providing quadrature coefficients via the trapezoidal method.
 */
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "quadrature_coeffs_nd.hpp"

/**
 * @brief Get the trapezoid coefficients in 1D.
 *
 * Calculate the quadrature coefficients for the trapezoid method defined on the provided index range.
 *
 * @param[in] idx_range
 * 	The index range on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the trapezoid method defined on the provided index range.
 */
template <class Grid>
host_t<FieldMem<double, IdxRange<Grid>>> trapezoid_quadrature_coefficients_1d(
        IdxRange<Grid> const& idx_range)
{
    host_t<FieldMem<double, IdxRange<Grid>>> coefficients(idx_range);
    IdxRange<Grid> middle_idx_range = idx_range.remove(IdxStep<Grid>(1), IdxStep<Grid>(1));

    coefficients(idx_range.front()) = 0.5 * distance_at_right(idx_range.front());
    ddc::for_each(middle_idx_range, [&](Idx<Grid> const idx) {
        coefficients(idx) = 0.5 * (distance_at_left(idx) + distance_at_right(idx));
    });
    coefficients(idx_range.back()) = 0.5 * distance_at_left(idx_range.back());

    if constexpr (Grid::continuous_dimension_type::PERIODIC) {
        coefficients(idx_range.front()) += 0.5 * distance_at_left(idx_range.back());
        coefficients(idx_range.back()) += 0.5 * distance_at_right(idx_range.front());
    }

    return coefficients;
}

/**
 * @brief Get the trapezoid coefficients in ND.
 *
 * Calculate the quadrature coefficients for the trapezoid method defined on the provided index range.
 *
 * @param[in] idx_range
 * 	The index range on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the trapezoid method defined on the provided index range.
 */
template <class... ODims>
host_t<FieldMem<double, IdxRange<ODims...>>> trapezoid_quadrature_coefficients(
        IdxRange<ODims...> const& idx_range)
{
    return quadrature_coeffs_nd(
            idx_range,
            (std::function<host_t<FieldMem<double, IdxRange<ODims>>>(IdxRange<ODims>)>(
                    trapezoid_quadrature_coefficients_1d<ODims>))...);
}
