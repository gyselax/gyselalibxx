// SPDX-License-Identifier: MIT

/**
 * @file spline_quadrature.hpp
 * File providing quadrature coefficients via a spline quadrature.
 */

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

#include <sll/matrix.hpp>

#include "ddc_aliases.hpp"



namespace {
template <class Grid1D>
using CoefficientFieldMem1D_host = host_t<DFieldMem<IdxRange<Grid1D>>>;
}


/**
 * @brief Get the spline quadrature coefficients.
 *
 * To integrate a function with a spline quadrature, we use:
 *
 * @f$ \int_a^b f(x)dx
 * \simeq \sum_{i = 0}^{N_{\text{basis}} -1 } c_i  \int_a^b b_{i,d}()x dx @f$,
 *
 * which rewritten gives
 *
 * @f$ \int_a^b f(x)dx
 * \simeq \sum_{i = 0}^{N_{\text{basis}} - 1} q_i f_i @f$,
 *
 * with
 *  - @f$\{ f_i\}_i @f$ the values of the function at the interpolation points;
 *  - @f$ q = \{ q_i\}_i @f$ the quadrature coefficients we compute thanks to
 *  @f$ q B^T = I_b @f$,
 *      - with @f$ B @f$ the matrix of B-splines @f$ B_{ij} = b_{j,d}(x_i)@f$,
 *      - and @f$ I_b = \int_a^b b_{i,d}(x)dx @f$ the integrated B-splines.
 *
 * More details are given in Emily Bourne's thesis
 * "Non-Uniform Numerical Schemes for the Modelling of Turbulence
 * in the 5D GYSELA Code". December 2022.
 *
 *
 * @param[in] idx_range
 *      The index range where the functions we want to integrate
 *      are defined.
 * @param[in] builder
 *      The spline builder describing the way in which splines would be constructed.
 *
 * @return A chunk with the quadrature coefficients @f$ q @f$.
 */
template <class Grid1D, class SplineBuilder>
host_t<DFieldMem<IdxRange<Grid1D>>> spline_quadrature_coefficients_1d(
        IdxRange<Grid1D> const& idx_range,
        SplineBuilder const& builder)
{
    static_assert(
            SplineBuilder::s_nbc_xmin == 0,
            "The spline quadrature requires a builder which can construct the coefficients using "
            "only the values at the interpolation points.");
    static_assert(
            SplineBuilder::s_nbc_xmax == 0,
            "The spline quadrature requires a builder which can construct the coefficients using "
            "only the values at the interpolation points.");
    // Since spline builder quadrature coeffs are not available on device, need host allocated builder.
    // See https://github.com/CExA-project/ddc/issues/598
    static_assert(
            (std::is_same_v<typename SplineBuilder::memory_space, Kokkos::HostSpace>),
            "SplineBuilder must be host allocated.");

    DFieldMem<IdxRange<Grid1D>, ddc::KokkosAllocator<double, typename SplineBuilder::memory_space>>
            quadrature_coefficients(builder.interpolation_domain());
    std::tie(std::ignore, quadrature_coefficients, std::ignore) = builder.quadrature_coefficients();
    return ddc::create_mirror_and_copy(quadrature_coefficients[idx_range]);
}



/**
 * @brief Get the spline quadrature coefficients in ND from N 1D quadrature coefficient.
 *
 * Calculate the quadrature coefficients for the spline quadrature method defined on the provided index range.
 *
 * @param[in] idx_range
 *      The index range on which the coefficients will be defined.
 * @param[in] builders
 *      The spline builder used for the quadrature coefficients in the different dimensions.
 *
 * @return The coefficients which define the spline quadrature method in ND.
 */
template <class ExecSpace, class... DDims, class... SplineBuilders>
DFieldMem<IdxRange<DDims...>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
spline_quadrature_coefficients(
        IdxRange<DDims...> const& idx_range,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::continuous_dimension_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D_host<DDims>...> current_dim_coeffs(
            spline_quadrature_coefficients_1d(ddc::select<DDims>(idx_range), builders)...);

    // Allocate ND coefficients
    DFieldMem<IdxRange<DDims...>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients(idx_range);
    auto coefficients_host = ddc::create_mirror(get_field(coefficients));
    // Serial loop is used due to nvcc bug concerning functions with variadic template arguments
    // (see https://github.com/kokkos/kokkos/pull/7059)
    ddc::for_each(idx_range, [&](Idx<DDims...> const idim) {
        // multiply the 1D coefficients by one another
        coefficients_host(idim)
                = (std::get<CoefficientFieldMem1D_host<DDims>>(current_dim_coeffs)(
                           ddc::select<DDims>(idim))
                   * ... * 1);
    });
    ddc::parallel_deepcopy(coefficients, coefficients_host);
    return coefficients;
}
