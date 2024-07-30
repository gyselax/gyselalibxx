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
using CoefficientChunk1D = host_t<FieldMem<double, IdxRange<Grid>>>;
}


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

    // Vector of integrals of B-splines
    host_t<FieldMem<double, IdxRange<bsplines_type>>> integral_bsplines(
            get_spline_idx_range(builder));
    ddc::discrete_space<bsplines_type>().integrals(get_field(integral_bsplines));

    // Solve matrix equation
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            integral_bsplines_mirror_with_additional_allocation(
                    "integral_bsplines_mirror_with_additional_allocation",
                    builder.get_interpolation_matrix().required_number_of_rhs_rows(),
                    1);
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            integral_bsplines_mirror = Kokkos::
                    subview(integral_bsplines_mirror_with_additional_allocation,
                            std::pair<std::size_t, std::size_t> {0, integral_bsplines.size()},
                            0);
    Kokkos::deep_copy(integral_bsplines_mirror, integral_bsplines.allocation_kokkos_view());
    builder.get_interpolation_matrix()
            .solve(integral_bsplines_mirror_with_additional_allocation, true);
    Kokkos::deep_copy(integral_bsplines.allocation_kokkos_view(), integral_bsplines_mirror);

    host_t<FieldMem<double, IdxRange<Grid>>> coefficients(idx_range);

    // Coefficients of quadrature in integral_bsplines (values which would always be multiplied
    // by f'(x)=0 are removed
    IdxRange<bsplines_type> slice
            = get_spline_idx_range(builder)
                      .remove(IdxStep<bsplines_type> {nbc_xmin}, IdxStep<bsplines_type> {nbc_xmax});

    Kokkos::deep_copy(
            coefficients.allocation_kokkos_view(),
            integral_bsplines[slice].allocation_kokkos_view());

    return coefficients;
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
template <class... DDims, class... SplineBuilders>
host_t<FieldMem<double, IdxRange<DDims...>>> neumann_spline_quadrature_coefficients(
        IdxRange<DDims...> const& idx_range,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::continuous_dimension_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D<DDims>...> current_dim_coeffs(
            neumann_spline_quadrature_coefficients_1d(ddc::select<DDims>(idx_range), builders)...);

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
