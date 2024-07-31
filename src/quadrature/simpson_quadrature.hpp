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
 * @param[in] idx_range
 * 	The index range on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the Simpson method defined on the provided index range.
 */
template <class ExecSpace, class Grid>
FieldMem<double, IdxRange<Grid>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
simpson_quadrature_coefficients_1d(IdxRange<Grid> const& idx_range)
{
    FieldMem<double, IdxRange<Grid>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients_alloc(idx_range);
    auto const coefficients = get_field(coefficients_alloc);
    double const dx_l = distance_at_left(idx_range.back());
    double const dx_r = distance_at_right(idx_range.front());
    Kokkos::parallel_for(
            "bounds",
            Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int i) {
                coefficients(idx_range.front()) = 1. / 3. * dx_r;
                coefficients(idx_range.back()) = 1. / 3. * dx_l;
                for (auto it = idx_range.begin() + 1; it < idx_range.end() - 1; it += 2) {
                    ddc::DiscreteElement<Grid> idx = *it;
                    coefficients(idx) = 2. / 3. * (distance_at_left(idx) + distance_at_right(idx));
                    idx += ddc::DiscreteVector<Grid>(1);
                    coefficients(idx) = 1. / 3. * (distance_at_left(idx) + distance_at_right(idx));
                }
                if (Grid::continuous_dimension_type::PERIODIC) {
                    coefficients(idx_range.front()) += 2. / 3. * dx_l;
                    coefficients(idx_range.back()) += 2. / 3. * dx_r;
                }
            });
    return coefficients_alloc;
}