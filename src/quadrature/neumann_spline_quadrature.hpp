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



namespace {
template <class IDim>
using CoefficientChunk1D = ddc::Chunk<double, ddc::DiscreteDomain<IDim>>;
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
 * @param[in] domain
 *      The domain on which the splines quadrature will be carried out.
 * @param[in] builder
 *      The spline builder used for the quadrature coefficients.
 *
 * @return The quadrature coefficients for the method defined on the provided domain.
 */
template <class IDim, class SplineBuilder>
ddc::Chunk<double, ddc::DiscreteDomain<IDim>> neumann_spline_quadrature_coefficients_1d(
        ddc::DiscreteDomain<IDim> const& domain,
        SplineBuilder const& builder)
{
    constexpr int nbe_xmin = SplineBuilder::s_nbe_xmin;
    constexpr int nbe_xmax = SplineBuilder::s_nbe_xmax;
    static_assert(
            SplineBuilder::s_bc_xmin == ddc::BoundCond::HERMITE,
            "The neumann spline quadrature requires a builder which uses Hermite boundary "
            "conditions.");
    static_assert(
            SplineBuilder::s_bc_xmax == ddc::BoundCond::HERMITE,
            "The neumann spline quadrature requires a builder which uses Hermite boundary "
            "conditions.");
    static_assert(
            nbe_xmin == 1,
            "The neumann spline quadrature requires a builder which uses the value of the "
            "derivative.");
    static_assert(
            nbe_xmax == 1,
            "The neumann spline quadrature requires a builder which uses the value of the "
            "derivative.");

    using bsplines_type = typename SplineBuilder::bsplines_type;

    assert(domain.size() == ddc::discrete_space<bsplines_type>().nbasis() - nbe_xmin - nbe_xmax);

    // Vector of integrals of B-splines
    ddc::Chunk<double, ddc::DiscreteDomain<bsplines_type>> integral_bsplines(
            builder.bsplines_domain());
    ddc::discrete_space<bsplines_type>().integrals(integral_bsplines.span_view());

    // Solve matrix equation
    builder.get_interpolation_matrix().solve_transpose_inplace(
            integral_bsplines.allocation_mdspan());

    ddc::Chunk<double, ddc::DiscreteDomain<IDim>> coefficients(domain);

    // Coefficients of quadrature in integral_bsplines (values which would always be multiplied
    // by f'(x)=0 are removed
    ddc::DiscreteDomain<bsplines_type> slice
            = builder.bsplines_domain()
                      .remove(ddc::DiscreteVector<bsplines_type> {nbe_xmin},
                              ddc::DiscreteVector<bsplines_type> {nbe_xmax});

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
 * @param[in] domain
 *      The domain on which the coefficients will be defined.
 * @param[in] builders
 *      The spline builder used for the quadrature coefficients in the different dimensions.
 *
 * @return The coefficients which define the spline quadrature method in ND.
 */
template <class... DDims, class... SplineBuilders>
ddc::Chunk<double, ddc::DiscreteDomain<DDims...>> neumann_spline_quadrature_coefficients(
        ddc::DiscreteDomain<DDims...> const& domain,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::bsplines_type::tag_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D<DDims>...> current_dim_coeffs(
            neumann_spline_quadrature_coefficients_1d(ddc::select<DDims>(domain), builders)...);

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
