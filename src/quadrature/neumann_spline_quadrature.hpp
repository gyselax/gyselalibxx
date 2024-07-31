// SPDX-License-Identifier: MIT

/**
 * @file neumann_spline_quadrature.hpp
 * File providing quadrature coefficients via a spline quadrature.
 */

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/matrix.hpp>

#include "ddc_aliases.hpp"



namespace {
template <class Grid>
using CoefficientFieldMem1D_h = host_t<FieldMem<double, IdxRange<Grid>>>;
template <class Grid>
using CoefficientField1D_h = host_t<Field<double, IdxRange<Grid>>>;

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
template <class Grid, class SplineBuilder>
host_t<FieldMem<double, IdxRange<Grid>>> neumann_spline_quadrature_coefficients_1d(
        IdxRange<Grid> const& idx_range,
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

    using bsplines_type = typename SplineBuilder::bsplines_type;

    assert(idx_range.size() == ddc::discrete_space<bsplines_type>().nbasis() - nbc_xmin - nbc_xmax);

    FieldMem<
            double,
            IdxRange<Grid>,
            ddc::KokkosAllocator<double, typename SplineBuilder::memory_space>>
            quadrature_coefficients(builder.interpolation_domain());
    // Even if derivatives coefficients on boundaries are eventually non-zero,
    // they are ignored for 0-flux Neumann boundary condition because
    // they would always be multiplied by f'(x)=0
    std::tie(std::ignore, quadrature_coefficients, std::ignore) = builder.quadrature_coefficients();
    return ddc::create_mirror_and_copy(quadrature_coefficients[idx_range]);
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
FieldMem<double, IdxRange<DDims...>, ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
neumann_spline_quadrature_coefficients(
        IdxRange<DDims...> const& idx_range,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::continuous_dimension_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientFieldMem1D_h<DDims>...> current_dim_coeffs_alloc(
            neumann_spline_quadrature_coefficients_1d(ddc::select<DDims>(idx_range), builders)...);
    std::tuple<CoefficientField1D_h<DDims>...> current_dim_coeffs(
            get_field(std::get<CoefficientFieldMem1D_h<DDims>>(current_dim_coeffs_alloc))...);
    // Allocate ND coefficients
    FieldMem<
            double,
            IdxRange<DDims...>,
            ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
            coefficients_alloc(idx_range);
    auto coefficients_alloc_host = ddc::create_mirror_view(get_field(coefficients_alloc));
    auto coefficients = get_field(coefficients_alloc_host);

    ddc::for_each(idx_range, [&](Idx<DDims...> const idim) {
        // multiply the 1D coefficients by one another
        coefficients(idim)
                = (std::get<CoefficientField1D_h<DDims>>(current_dim_coeffs)(
                           ddc::select<DDims>(idim))
                   * ... * 1);
    });
    ddc::parallel_deepcopy(coefficients_alloc, coefficients_alloc_host);
    return coefficients_alloc;
}
