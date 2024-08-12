// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

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
template <class ExecSpace, class Grid1D>
FieldMem<double, IdxRange<Grid1D>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
simpson_quadrature_coefficients_1d(IdxRange<Grid1D> const& idx_range)
{
    DFieldMem<IdxRange<Grid1D>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients_alloc(idx_range);
    DField<IdxRange<Grid1D>,
           std::experimental::layout_right,
           typename ExecSpace::memory_space> const coefficients
            = get_field(coefficients_alloc);
    double const dx_l = distance_at_left(idx_range.back());
    double const dx_r = distance_at_right(idx_range.front());

    Kokkos::parallel_for(
            "bounds",
            Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int i) {
                coefficients(idx_range.front()) = 1. / 3. * dx_r;
                coefficients(idx_range.back()) = 1. / 3. * dx_l;
            });
    Kokkos::parallel_for(
            "centre",
            Kokkos::RangePolicy<
                    ExecSpace>(0, (*idx_range.end() - *idx_range.begin()).value() / 2 - 1),
            KOKKOS_LAMBDA(const int i) {
                Idx<Grid1D> idx = *idx_range.begin() + 1 + 2 * i;
                coefficients(idx) = 2. / 3. * (distance_at_left(idx) + distance_at_right(idx));
                idx += IdxStep<Grid1D>(1);
                coefficients(idx) = 1. / 3. * (distance_at_left(idx) + distance_at_right(idx));
            });
    if constexpr (Grid1D::continuous_dimension_type::PERIODIC) {
        Kokkos::parallel_for(
                "bounds",
                Kokkos::RangePolicy<ExecSpace>(0, 1),
                KOKKOS_LAMBDA(const int i) {
                    coefficients(idx_range.front()) += 2. / 3. * dx_l;
                    coefficients(idx_range.back()) += 2. / 3. * dx_r;
                });
    }


    return coefficients_alloc;
}