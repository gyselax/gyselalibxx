// SPDX-License-Identifier: MIT
/**
 * @file quadrature_coeffs_nd.hpp
 * File providing helper functions for defining multi-dimensional quadrature methods.
 */
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


namespace {
template <class ExecSpace, class Grid1D>
using CoefficientFieldMem1D = DFieldMem<
        IdxRange<Grid1D>,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>;
template <class ExecSpace, class Grid1D>
using CoefficientField1D = DField<
        IdxRange<Grid1D>,
        std::experimental::layout_right,
        typename ExecSpace::memory_space>;
} // namespace

/**
 * @brief Helper function which creates ND dimensions from N 1D quadrature coefficient functions.
 *
 * @param idx_range
 *      The index range on which the coefficients will be defined.
 * @param funcs
 *      The functions which define quadrature coefficients in the different dimensions.
 *
 * @returns The coefficients which define the quadrature method in ND.
 */
template <class ExecSpace, class... DDims>
DFieldMem<IdxRange<DDims...>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
quadrature_coeffs_nd(
        IdxRange<DDims...> const& idx_range,
        std::function<DFieldMem<
                IdxRange<DDims>,
                ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>(
                IdxRange<DDims>)>... funcs)
{
    DFieldMem<IdxRange<DDims...>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients_alloc(idx_range);
    DField<IdxRange<DDims...>, std::experimental::layout_right, typename ExecSpace::memory_space>
            coefficients(get_field(coefficients_alloc));
    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D<ExecSpace, DDims>...> current_dim_coeffs_alloc(
            funcs(ddc::select<DDims>(idx_range))...);
    std::tuple<CoefficientField1D<ExecSpace, DDims>...> current_dim_coeffs(get_field(
            std::get<CoefficientFieldMem1D<ExecSpace, DDims>>(current_dim_coeffs_alloc))...);

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range,
            KOKKOS_LAMBDA(Idx<DDims...> const idim) {
                // multiply the 1D coefficients by one another

                coefficients(idim)
                        = (std::get<CoefficientField1D<ExecSpace, DDims>>(current_dim_coeffs)(
                                   ddc::select<DDims>(idim))
                           * ... * 1);
            });
    return coefficients_alloc;
}
