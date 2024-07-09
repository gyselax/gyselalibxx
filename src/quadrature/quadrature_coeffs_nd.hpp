// SPDX-License-Identifier: MIT
/**
 * @file quadrature_coeffs_nd.hpp
 * File providing helper functions for defining multi-dimensional quadrature methods.
 */
#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>


namespace {
template <class ExecSpace, class IDim>
using CoefficientChunk1D = ddc::Chunk<
        double,
        ddc::DiscreteDomain<IDim>,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>;
template <class ExecSpace, class IDim>
using CoefficientChunkSpan1D = ddc::ChunkSpan<
        double,
        ddc::DiscreteDomain<IDim>,
        std::experimental::layout_right,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>;
} // namespace

/**
 * @brief Helper function which creates ND dimensions from N 1D quadrature coefficient functions.
 *
 * @param domain
 *      The domain on which the coefficients will be defined.
 * @param funcs
 *      The functions which define quadrature coefficients in the different dimensions.
 *
 * @returns The coefficients which define the quadrature method in ND.
 */
template <class ExecSpace, class... DDims>
ddc::Chunk<
        double,
        ddc::DiscreteDomain<DDims...>,
        ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
quadrature_coeffs_nd(
        ddc::DiscreteDomain<DDims...> const& domain,
        std::function<ddc::Chunk<
                double,
                ddc::DiscreteDomain<DDims>,
                ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>(
                ddc::DiscreteDomain<DDims>)>... funcs)
{
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<DDims...>>> coefficients_alloc(domain);
    ddc::ChunkSpan coefficients = coefficients_alloc.span_view();
    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D<ExecSpace, DDims>...> current_dim_coeffs_alloc(
            funcs(ddc::select<DDims>(domain))...);
    std::tuple<CoefficientChunkSpan1D<ExecSpace, DDims>...> current_dim_coeffs(
            std::get<CoefficientChunk1D<ExecSpace, DDims>>(current_dim_coeffs_alloc)
                    .span_view()...);

    ddc::parallel_for_each(
            ExecSpace(),
            domain,
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDims...> const idim) {
                // multiply the 1D coefficients by one another

                coefficients(idim)
                        = (std::get<CoefficientChunkSpan1D<ExecSpace, DDims>>(current_dim_coeffs)(
                                   ddc::select<DDims>(idim))
                           * ... * 1);
            });
    return std::move(coefficients_alloc);
}
