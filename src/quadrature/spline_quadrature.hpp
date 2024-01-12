// SPDX-License-Identifier: MIT

/**
 * @file spline_quadrature.hpp
 * File providing quadrature coefficients via a spline quadrature.
 */

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/matrix.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>



namespace {
template <class IDim>
using CoefficientChunk1D = ddc::Chunk<double, ddc::DiscreteDomain<IDim>>;
}


/**
 * @brief Get the spline quadrature coefficients in 1D.
 *
 * This function mainly uses the SplineBuilder::quadrature_coefficients function.
 *
 * @param[in] domain
 *      The domain on which the splines quadrature will be carried out.
 * @param[in] builder
 *      The spline builder used for the quadrature coefficients.
 *
 * @return The quadrature coefficients for the method defined on the provided domain.
 */
template <class IDim, class SplineBuilder>
ddc::Chunk<double, ddc::DiscreteDomain<IDim>> spline_quadrature_coefficients_1d(
        ddc::DiscreteDomain<IDim> const& domain,
        SplineBuilder const& builder)
{
    static_assert(
            SplineBuilder::s_nbe_xmin == 0,
            "The spline quadrature requires a builder which can construct the coefficients using "
            "only the values at the interpolation points.");
    static_assert(
            SplineBuilder::s_nbe_xmax == 0,
            "The spline quadrature requires a builder which can construct the coefficients using "
            "only the values at the interpolation points.");
    return builder.quadrature_coefficients(domain);
}



/**
 * @brief Get the spline quadrature coefficients in ND from N 1D quadrature coefficient.
 *
 * Calculate the quadrature coefficients for the spline quadrature method defined on the provided domain.
 *
 * @param[in] domain
 *      The domain on which the coefficients will be defined.
 * @param[in] builders
 *      The spline builder used for the quadrature coefficients in the different dimensions.
 *
 * @return The coefficients which define the spline quadrature method in ND.
 */
template <class... DDims, class... SplineBuilders>
ddc::Chunk<double, ddc::DiscreteDomain<DDims...>> spline_quadrature_coefficients(
        ddc::DiscreteDomain<DDims...> const& domain,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::bsplines_type::tag_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D<DDims>...> current_dim_coeffs(
            spline_quadrature_coefficients_1d(ddc::select<DDims>(domain), builders)...);

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
