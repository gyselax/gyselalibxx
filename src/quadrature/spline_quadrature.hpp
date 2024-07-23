// SPDX-License-Identifier: MIT

/**
 * @file spline_quadrature.hpp
 * File providing quadrature coefficients via a spline quadrature.
 */

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

#include <sll/matrix.hpp>



namespace {
template <class IDim>
using CoefficientChunk1D_h = ddc::Chunk<double, ddc::DiscreteDomain<IDim>>;
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
 * @param[in] domain
 *      The domain where the functions we want to integrate
 *      are defined.
 * @param[in] builder
 *      The spline builder describing the way in which splines would be constructed.
 *
 * @return A chunk with the quadrature coefficients @f$ q @f$.
 */
template <class IDim, class SplineBuilder>
ddc::Chunk<double, ddc::DiscreteDomain<IDim>> spline_quadrature_coefficients_1d(
        ddc::DiscreteDomain<IDim> const& domain,
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

    using bsplines_type = typename SplineBuilder::bsplines_type;

    // Vector of integrals of B-splines
    ddc::Chunk<double, ddc::DiscreteDomain<bsplines_type>> integral_bsplines(
            builder.spline_domain());
    ddc::discrete_space<bsplines_type>().integrals(integral_bsplines.span_view());

    // Solve matrix equation
    ddc::ChunkSpan integral_bsplines_without_periodic_point
            = integral_bsplines.span_view()[ddc::DiscreteDomain<bsplines_type>(
                    ddc::DiscreteElement<bsplines_type>(0),
                    ddc::DiscreteVector<bsplines_type>(builder.get_interpolation_matrix().size()))];
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            integral_bsplines_mirror_with_additional_allocation(
                    "integral_bsplines_mirror_with_additional_allocation",
                    builder.get_interpolation_matrix().required_number_of_rhs_rows(),
                    1);
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            integral_bsplines_mirror = Kokkos::
                    subview(integral_bsplines_mirror_with_additional_allocation,
                            std::pair<std::size_t, std::size_t> {
                                    0,
                                    integral_bsplines_without_periodic_point.size()},
                            0);
    Kokkos::deep_copy(
            integral_bsplines_mirror,
            integral_bsplines_without_periodic_point.allocation_kokkos_view());

    // Solve matrix equation
    builder.get_interpolation_matrix()
            .solve(integral_bsplines_mirror_with_additional_allocation, true);
    Kokkos::deep_copy(
            integral_bsplines_without_periodic_point.allocation_kokkos_view(),
            integral_bsplines_mirror);

    ddc::Chunk<double, ddc::DiscreteDomain<IDim>> coefficients(domain);
    Kokkos::deep_copy(
            coefficients.allocation_kokkos_view(),
            integral_bsplines_without_periodic_point.allocation_kokkos_view());

    return coefficients;
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
device_t<ddc::Chunk<double, ddc::DiscreteDomain<DDims...>>> spline_quadrature_coefficients(
        ddc::DiscreteDomain<DDims...> const& domain,
        SplineBuilders const&... builders)
{
    assert((std::is_same_v<
                    typename DDims::continuous_dimension_type,
                    typename SplineBuilders::continuous_dimension_type> and ...));

    // Get coefficients for each dimension
    std::tuple<CoefficientChunk1D_h<DDims>...> current_dim_coeffs(
            spline_quadrature_coefficients_1d(ddc::select<DDims>(domain), builders)...);

    // Allocate ND coefficients
    ddc::Chunk<double, ddc::DiscreteDomain<DDims...>> coefficients_host(domain);
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<DDims...>>> coefficients(domain);

    ddc::for_each(domain, [&](ddc::DiscreteElement<DDims...> const idim) {
        // multiply the 1D coefficients by one another
        coefficients_host(idim)
                = (std::get<CoefficientChunk1D_h<DDims>>(current_dim_coeffs)(
                           ddc::select<DDims>(idim))
                   * ... * 1);
    });
    ddc::parallel_deepcopy(coefficients, coefficients_host);
    return coefficients;
}
