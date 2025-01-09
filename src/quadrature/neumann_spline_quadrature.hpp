// SPDX-License-Identifier: MIT

/**
 * @file neumann_spline_quadrature.hpp
 * File providing quadrature coefficients via a spline quadrature.
 */

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"



namespace {
template <class ExecSpace, class Grid1D>
using CoefficientFieldMem1D = DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space>;
template <class ExecSpace, class Grid1D>
using CoefficientField1D = DField<IdxRange<Grid1D>, typename ExecSpace::memory_space>;

} // namespace



/**
 * @brief Get the spline quadrature coefficients in 1D.
 *
 * This function calculates the quadrature coefficients which define a quadrature equivalent
 * to calculating and integrating a spline approximation of a function. The spline approximation
 * would be calculated with homogeneous Neumann boundary conditions.
 * This method of defining quadrature coefficients is described in section Emily Bourne's thesis[1].
 *
 * [1] Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GYSELA Code
 *     Emily Bourne, December 2022-
 *
 * @param[in] idx_range
 *      The index range on which the splines quadrature will be carried out.
 * @param[in] builder
 *      The spline builder used for the quadrature coefficients.
 *
 * @return The quadrature coefficients for the method defined on the provided index range.
 */
template <class ExecSpace, class Grid1D, class SplineBuilder>
DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space>
neumann_spline_quadrature_coefficients_1d(
        IdxRange<Grid1D> const& idx_range,
        SplineBuilder const& builder)
{
    constexpr int nbc_xmin = SplineBuilder::s_nbc_xmin;
    constexpr int nbc_xmax = SplineBuilder::s_nbc_xmax;
    static_assert(
            SplineBuilder::s_bc_xmin == ddc::BoundCond::HERMITE,
            "The neumann spline quadrature requires a builder which uses Hermite boundary "
            "conditions.");
    static_assert(
            SplineBuilder::s_bc_xmax == ddc::BoundCond::HERMITE,
            "The neumann spline quadrature requires a builder which uses Hermite boundary "
            "conditions.");
    static_assert(
            nbc_xmin == 1,
            "The neumann spline quadrature requires a builder which uses the value of the "
            "derivative.");
    static_assert(
            nbc_xmax == 1,
            "The neumann spline quadrature requires a builder which uses the value of the "
            "derivative.");
    assert(idx_range.size()
           == ddc::discrete_space<typename SplineBuilder::bsplines_type>().nbasis() - nbc_xmin
                      - nbc_xmax);

    DFieldMem<IdxRange<Grid1D>, typename SplineBuilder::memory_space> quadrature_coefficients(
            builder.interpolation_domain());
    // Even if derivatives coefficients on boundaries are eventually non-zero,
    // they are ignored for 0-flux Neumann boundary condition because
    // they would always be multiplied by f'(x)=0
    std::tie(std::ignore, quadrature_coefficients, std::ignore) = builder.quadrature_coefficients();
    DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space> output_quad_coefficients(
            idx_range);
    ddc::parallel_deepcopy(output_quad_coefficients, quadrature_coefficients[idx_range]);
    return output_quad_coefficients;
}



/**
 * @brief Get the spline quadrature coefficients in ND from N 1D quadrature coefficient.
 *
 * This function calculates the quadrature coefficients which define a quadrature equivalent
 * to calculating and integrating a spline approximation of a function. The spline approximation
 * would be calculated with homogeneous Neumann boundary conditions.
 *
 * @param[in] idx_range
 *      The index range on which the coefficients will be defined.
 * @param[in] builders
 *      The spline builder used for the quadrature coefficients in the different dimensions.
 *
 * @return The coefficients which define the spline quadrature method in ND.
 */
template <class ExecSpace, class... DDims, class... SplineBuilders>
DFieldMem<IdxRange<DDims...>, typename ExecSpace::memory_space>
neumann_spline_quadrature_coefficients(
        IdxRange<DDims...> const& idx_range,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::continuous_dimension_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D<Kokkos::DefaultHostExecutionSpace, DDims>...>
    current_dim_coeffs_alloc(
            neumann_spline_quadrature_coefficients_1d<
                    Kokkos::DefaultHostExecutionSpace>(ddc::select<DDims>(idx_range), builders)...);
    std::tuple<CoefficientField1D<Kokkos::DefaultHostExecutionSpace, DDims>...> current_dim_coeffs(
            get_field(std::get<CoefficientFieldMem1D<Kokkos::DefaultHostExecutionSpace, DDims>>(
                    current_dim_coeffs_alloc))...);
    // Allocate ND coefficients
    DFieldMem<IdxRange<DDims...>, typename ExecSpace::memory_space> coefficients_alloc(idx_range);
    auto coefficients_alloc_host = ddc::create_mirror(get_field(coefficients_alloc));
    host_t<DField<IdxRange<DDims...>>> coefficients(get_field(coefficients_alloc_host));
    // Serial loop is used due to nvcc bug concerning functions with variadic template arguments
    // (see https://github.com/kokkos/kokkos/pull/7059)
    ddc::for_each(idx_range, [&](Idx<DDims...> const idim) {
        // multiply the 1D coefficients by one another
        coefficients(idim)
                = (std::get<CoefficientField1D<Kokkos::DefaultHostExecutionSpace, DDims>>(
                           current_dim_coeffs)(ddc::select<DDims>(idim))
                   * ... * 1);
    });
    ddc::parallel_deepcopy(coefficients_alloc, coefficients_alloc_host);
    return coefficients_alloc;
}
