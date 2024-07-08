// SPDX-License-Identifier: MIT
/**
 * @file quadrature_coeffs_nd.hpp
 * File providing helper functions for defining multi-dimensional quadrature methods.
 */
#pragma once

#include <ddc/ddc.hpp>


namespace {
template <class IDim>
using CoefficientChunk1D = device_t<ddc::Chunk<double, ddc::DiscreteDomain<IDim>>>;
template <class IDim>
using CoefficientChunkSpan1D = device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<IDim>>>;
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
device_t<ddc::Chunk<double, ddc::DiscreteDomain<DDims...>>> quadrature_coeffs_nd(
        ddc::DiscreteDomain<DDims...> const& domain,
        std::function<device_t<ddc::Chunk<double, ddc::DiscreteDomain<DDims>>>(
                ddc::DiscreteDomain<DDims>)>... funcs)
{
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<DDims...>>> coefficients_alloc(domain);
    ddc::ChunkSpan coefficients = coefficients_alloc.span_view();
    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D<DDims>...> current_dim_coeffs_alloc(
            funcs(ddc::select<DDims>(domain))...);
    std::tuple<CoefficientChunkSpan1D<DDims>...> current_dim_coeffs(
            std::get<CoefficientChunk1D<DDims>>(current_dim_coeffs_alloc).span_view()...);

    ddc::parallel_for_each(
            ExecSpace(),
            domain,
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDims...> const idim) {
                // multiply the 1D coefficients by one another

                coefficients(idim)
                        = (std::get<CoefficientChunkSpan1D<DDims>>(current_dim_coeffs)(
                                   ddc::select<DDims>(idim))
                           * ... * 1);
            });
    return std::move(coefficients_alloc);
}
