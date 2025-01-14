// SPDX-License-Identifier: MIT
/**
 * @file volume_quadrature_nd.hpp
 * File providing functions to adapt quadrature coefficients so they can integrate over a N-D volume.
 * A N-D volume is a surface in 2D and a volume in 3D.
 */

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "quadrature.hpp"



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
 * @param[in] exec_space
 *     The space on which the function is executed (CPU/GPU).
 * @param[in] mapping
 *      The mapping function from the logical index range @f$ (r,\theta) @f$
 *      to the physical index range @f$ (x, y) @f$.
 * @param[in] coefficients_alloc
 *      The quadrature coefficients @f$\{q_{ij}\}_{ij} @f$.
 *
 * @return A rvalue FieldMem to the modified coefficients  @f$\{q_{ij}| det(J(r_i,\theta_j))|\}_{ij} @f$.
 */
template <class Mapping, class IdxRangeCoeffs, class ExecSpace>
DFieldMem<IdxRangeCoeffs, typename ExecSpace::memory_space> compute_coeffs_on_mapping(
        ExecSpace exec_space,
        Mapping& mapping,
        DFieldMem<IdxRangeCoeffs, typename ExecSpace::memory_space>&& coefficients_alloc)
{
    static_assert(is_curvilinear_2d_mapping_v<Mapping>);

    using R = typename Mapping::curvilinear_tag_r;
    using Theta = typename Mapping::curvilinear_tag_theta;

    static_assert(has_2d_jacobian_v<Mapping, Coord<R, Theta>>);

    using IdxCoeffs = typename IdxRangeCoeffs::discrete_element_type;
    IdxRangeCoeffs grid = get_idx_range(coefficients_alloc);
    DField<IdxRangeCoeffs, typename ExecSpace::memory_space> coefficients(coefficients_alloc);
    ddc::parallel_for_each(
            exec_space,
            grid,
            KOKKOS_LAMBDA(IdxCoeffs const idx) {
                Coord<R, Theta> coord(ddc::coordinate(idx));
                coefficients(idx) *= fabs(mapping.jacobian(coord));
            });
    return std::move(coefficients_alloc);
}
