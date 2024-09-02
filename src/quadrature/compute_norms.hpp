// SPDX-License-Identifier: MIT
/**
 * @file compute_norms.hpp
 * File providing the L1 and the L2 norms.
 */

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "quadrature.hpp"


/**
 * @brief Compute L1 norm of a function with a given quadrature.
 *
 * @f$ \int_{\Omega} |f(X)|  dX @f$
 *
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      A Field to the value of the function on the quadrature grid.
 *
 * @return A double containing the L1 norm of the function.
 */
template <class IdxRange>
double compute_L1_norm(
        host_t<Quadrature<IdxRange>> quadrature,
        Field<double,
              IdxRange,
              std::experimental::layout_right,
              Kokkos::DefaultHostExecutionSpace::memory_space> function)
{
    using Idx = typename IdxRange::discrete_element_type;
    return quadrature(
            Kokkos::DefaultHostExecutionSpace(),
            KOKKOS_LAMBDA(Idx const idx) { return Kokkos::fabs(function(idx)); });
}



/**
 * @brief Compute L2 norm of a function with a given quadrature.
 *
 * @f$ \sqrt{\int_{\Omega} |f(X)|^2  dX} @f$
 *
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      A Field to the value of the function on the quadrature grid.
 *
 * @return A double containing the L2 norm of the function.
 */
template <class IdxRange>
double compute_L2_norm(
        host_t<Quadrature<IdxRange>> quadrature,
        Field<double,
              IdxRange,
              std::experimental::layout_right,
              Kokkos::DefaultHostExecutionSpace::memory_space> function)
{
    using Idx = typename IdxRange::discrete_element_type;
    return std::sqrt(quadrature(
            Kokkos::DefaultHostExecutionSpace(),
            KOKKOS_LAMBDA(Idx const idx) { return function(idx) * function(idx); }));
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
 *      The mapping function from the logical index range @f$ (r,\theta) @f$
 *      to the physical index range @f$ (x, y) @f$.
 * @param[in] coefficients
 *      The quadrature coefficients @f$\{q_{ij}\}_{ij} @f$.
 *
 * @return A rvalue FieldMem to the modified coefficients  @f$\{q_{ij}| det(J(r_i,\theta_j))|\}_{ij} @f$.
 */
template <class Mapping, class... IDim>
host_t<FieldMem<double, IdxRange<IDim...>>> compute_coeffs_on_mapping(
        Mapping& mapping,
        FieldMem<
                double,
                IdxRange<IDim...>,
                ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>&&
                coefficients)
{
    IdxRange<IDim...> grid = get_idx_range<IDim...>(coefficients);
    ddc::for_each(grid, [&](Idx<IDim...> const idx) {
        Coord<typename Mapping::curvilinear_tag_r, typename Mapping::curvilinear_tag_theta> coord(
                ddc::coordinate(idx));
        coefficients(idx) *= fabs(mapping.jacobian(coord));
    });
    return std::move(coefficients);
}
