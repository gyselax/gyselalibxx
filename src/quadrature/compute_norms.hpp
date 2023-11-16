/**
 * @file compute_norms.hpp
 * File providing the L1 and the L2 norms.
 */

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "quadrature.hpp"


/**
 * @brief Compute L1 norm of a function with a given quadrature.
 *
 * @f$ \int_{\Omega} |f(X)|  dX @f$
 *
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      A ChunkSpan to the value of the function on the quadrature grid.
 *
 * @return A double containing the L1 norm of the function.
 */
template <class... IDim>
double compute_L1_norm(
        Quadrature<IDim...> quadrature,
        ddc::ChunkSpan<double, ddc::DiscreteDomain<IDim...>> function)
{
    auto grid = ddc::get_domain<IDim...>(function);
    ddc::Chunk<double, ddc::DiscreteDomain<IDim...>> abs_function(grid);
    ddc::for_each(grid, [&](ddc::DiscreteElement<IDim...> const idx) {
        abs_function(idx) = fabs(function(idx));
    });

    return quadrature(abs_function);
}



/**
 * @brief Compute L2 norm of a function with a given quadrature.
 *
 * @f$ \sqrt{\int_{\Omega} |f(X)|^2  dX} @f$
 *
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      A ChunkSpan to the value of the function on the quadrature grid.
 *
 * @return A double containing the L2 norm of the function.
 */
template <class... IDim>
double compute_L2_norm(
        Quadrature<IDim...> quadrature,
        ddc::ChunkSpan<double, ddc::DiscreteDomain<IDim...>> function)
{
    auto grid = ddc::get_domain<IDim...>(function);
    ddc::Chunk<double, ddc::DiscreteDomain<IDim...>> square_function(grid);
    ddc::for_each(grid, [&](ddc::DiscreteElement<IDim...> const idx) {
        square_function(idx) = function(idx) * function(idx);
    });

    return std::sqrt(quadrature(square_function));
}



/**
 * @brief Add the Jacobian determinant to the coefficients.
 *
 * For polar integration, we can add the Jacobian determinant to the
 * quadrature coefficients.
 *
 * - 2D example: (but it is implemented for ND)
 *
 * @f$ \int_{\Omega_{\text{cart}}} f(x,y) dxdy
 *  = \int_{\Omega} f(r,\theta) |det(J(r,\theta))| drd\theta
 *  \simeq \sum_{i,j} [q_{i,j}| det(J(r_i,\theta_j))|] f_{ij}@f$
 *
 * This function uses rvalues. It means that coefficients is a temporary input
 * parameter and it returns a temporary coefficient object. The Quadrature
 * object can only be instantiate with rvalues.
 *
 * @param[in] mapping
 *      The mapping function from the logical domain @f$ (r,\theta) @f$
 *      to the physical domain @f$ (x, y) @f$.
 * @param[in, out] coefficients
 *      The quadrature coefficients @f$\{q_{ij}\}_{ij} @f$.
 *
 * @return A rvalue Chunk to the modified coefficients  @f$\{q_{ij}| det(J(r_i,\theta_j))|\}_{ij} @f$.
 */
template <class Mapping, class... IDim>
ddc::Chunk<double, ddc::DiscreteDomain<IDim...>> compute_coeffs_on_mapping(
        Mapping& mapping,
        ddc::Chunk<double, ddc::DiscreteDomain<IDim...>>&& coefficients)
{
    ddc::DiscreteDomain<IDim...> grid = ddc::get_domain<IDim...>(coefficients);
    ddc::for_each(grid, [&](ddc::DiscreteElement<IDim...> const idx) {
        ddc::Coordinate<typename Mapping::curvilinear_tag_r, typename Mapping::curvilinear_tag_p>
                coord(ddc::coordinate(idx));
        coefficients(idx) *= fabs(mapping.jacobian(coord));
    });
    return std::move(coefficients);
}
