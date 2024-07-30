// SPDX-License-Identifier: MIT
/**
 * @file quadrature_coeffs_nd.hpp
 * File providing helper functions for defining multi-dimensional quadrature methods.
 */
#pragma once

#include <ddc/ddc.hpp>


namespace {
template <class IDim>
using CoefficientChunk1D = ddc::Chunk<double, ddc::DiscreteDomain<IDim>>;
}

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
template <class... DDims>
ddc::Chunk<double, ddc::DiscreteDomain<DDims...>> quadrature_coeffs_nd(
        ddc::DiscreteDomain<DDims...> const& domain,
        std::function<ddc::Chunk<double, ddc::DiscreteDomain<DDims>>(
                ddc::DiscreteDomain<DDims>)>... funcs)
{
    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D<DDims>...> current_dim_coeffs(
            funcs(ddc::select<DDims>(domain))...);

    // Allocate ND coefficients
    ddc::Chunk<double, ddc::DiscreteDomain<DDims...>> coefficients(domain);

    ddc::for_each(domain, [&](ddc::DiscreteElement<DDims...> const idim) {
        // multiply the 1D coefficients by one another
        coefficients(idim)
                = (std::get<CoefficientChunk1D<DDims>>(current_dim_coeffs)(ddc::select<DDims>(idim))
                   * ... * 1);
    });

    return coefficients;
}
