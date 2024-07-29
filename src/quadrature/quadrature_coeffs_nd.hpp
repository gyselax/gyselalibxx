// SPDX-License-Identifier: MIT
/**
 * @file quadrature_coeffs_nd.hpp
 * File providing helper functions for defining multi-dimensional quadrature methods.
 */
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


namespace {
template <class Grid>
using CoefficientChunk1D = host_t<FieldMem<double, IdxRange<Grid>>>;
}

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
template <class... DDims>
host_t<FieldMem<double, IdxRange<DDims...>>> quadrature_coeffs_nd(
        IdxRange<DDims...> const& idx_range,
        std::function<host_t<FieldMem<double, IdxRange<DDims>>>(IdxRange<DDims>)>... funcs)
{
    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D<DDims>...> current_dim_coeffs(
            funcs(ddc::select<DDims>(idx_range))...);

    // Allocate ND coefficients
    host_t<FieldMem<double, IdxRange<DDims...>>> coefficients(idx_range);

    ddc::for_each(idx_range, [&](Idx<DDims...> const idim) {
        // multiply the 1D coefficients by one another
        coefficients(idim)
                = (std::get<CoefficientChunk1D<DDims>>(current_dim_coeffs)(ddc::select<DDims>(idim))
                   * ... * 1);
    });

    return coefficients;
}
