// SPDX-License-Identifier: MIT
/**
 * @file quadrature_coeffs_nd.hpp
 * File providing helper functions for defining multi-dimensional quadrature methods.
 */
#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>


namespace {
template <class ExecSpace, class Grid>
using CoefficientFieldMem1D = FieldMem<
        double,
        IdxRange<Grid>,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>;
template <class ExecSpace, class Grid>
using CoefficientField1D
        = Field<double,
                IdxRange<Grid>,
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
FieldMem<double, IdxRange<DDims...>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
quadrature_coeffs_nd(
        IdxRange<DDims...> const& idx_range,
        std::function<FieldMem<
                double,
                IdxRange<DDims>,
                ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>(
                IdxRange<DDims>)>... funcs)
{
    ddc::Chunk<
            double,
            ddc::DiscreteDomain<DDims...>,
            ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients_alloc(idx_range);
    ddc::ChunkSpan coefficients = get_field(coefficients_alloc);
    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D<ExecSpace, DDims>...> current_dim_coeffs_alloc(
            funcs(ddc::select<DDims>(idx_range))...);
    std::tuple<CoefficientField1D<ExecSpace, DDims>...> current_dim_coeffs(get_field(
            std::get<CoefficientFieldMem1D<ExecSpace, DDims>>(current_dim_coeffs_alloc))...);

    ddc::parallel_for_each(
            ExecSpace(),
            idx_range,
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDims...> const idim) {
                // multiply the 1D coefficients by one another

                coefficients(idim)
                        = (std::get<CoefficientField1D<ExecSpace, DDims>>(current_dim_coeffs)(
                                   ddc::select<DDims>(idim))
                           * ... * 1);
            });
    return coefficients_alloc;
}
