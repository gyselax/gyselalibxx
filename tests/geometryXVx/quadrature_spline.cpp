// SPDX-License-Identifier: MIT

#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/spline_boundary_conditions.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.hpp"
#include "quadrature.hpp"
#include "spline_quadrature.hpp"

namespace {

TEST(SplineQuadratureTest, ExactForConstantFunc)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IVectX const x_size(10);

    // Creating mesh & supports
    using sllBSplinesX = UniformBSplines<RDimX, 3>;
    auto constexpr sllSplineXBoundary = RDimX::PERIODIC ? BoundCond::PERIODIC : BoundCond::GREVILLE;
    using sllGrevillePointsX
            = GrevilleInterpolationPoints<sllBSplinesX, sllSplineXBoundary, sllSplineXBoundary>;
    using sllIDimX = typename sllGrevillePointsX::interpolation_mesh_type;
    using sllSplineXBuilder
            = SplineBuilder<sllBSplinesX, sllIDimX, sllSplineXBoundary, sllSplineXBoundary>;
    using sllIDomainX = ddc::DiscreteDomain<sllIDimX>;
    using sllDFieldX = device_t<ddc::Chunk<double, sllIDomainX>>;

    ddc::init_discrete_space<sllBSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<sllIDimX>(sllGrevillePointsX::get_sampling());
    sllIDomainX interpolation_domain_x(sllGrevillePointsX::get_domain());

    sllSplineXBuilder const builder_x(interpolation_domain_x);

    sllIDomainX const gridx = builder_x.interpolation_domain();

    host_t<sllDFieldX> const quadrature_coeffs = spline_quadrature_coefficients(gridx, builder_x);
    Quadrature const integrate(quadrature_coeffs.span_cview());

    host_t<sllDFieldX> values(gridx);

    ddc::for_each(gridx, [&](ddc::DiscreteElement<IDimX> const idx) { values(idx) = 1.0; });
    double integral = integrate(values);
    double expected_val = x_max - x_min;
    EXPECT_LE(abs(integral - expected_val), 1e-15);
}


template <std::size_t N>
struct Y
{
    static bool constexpr PERIODIC = false;
};

template <std::size_t N>
double compute_error(int n_elems)
{
    using DimY = Y<N>;
    using BSplinesY = UniformBSplines<DimY, 3>;
    using GrevillePointsY
            = GrevilleInterpolationPoints<BSplinesY, BoundCond::GREVILLE, BoundCond::GREVILLE>;
    using IDimY = typename GrevillePointsY::interpolation_mesh_type;
    using SplineYBuilder
            = SplineBuilder<BSplinesY, IDimY, BoundCond::GREVILLE, BoundCond::GREVILLE>;
    using IDomainY = ddc::DiscreteDomain<IDimY>;
    using DFieldY = device_t<ddc::Chunk<double, IDomainY>>;

    ddc::Coordinate<Y<N>> const y_min(0.0);
    ddc::Coordinate<Y<N>> const y_max(M_PI);

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, n_elems);

    ddc::init_discrete_space<IDimY>(GrevillePointsY::get_sampling());
    IDomainY const gridy(GrevillePointsY::get_domain());

    SplineYBuilder const builder_y(gridy);

    host_t<DFieldY> const quadrature_coeffs = spline_quadrature_coefficients(gridy, builder_y);
    Quadrature const integrate(quadrature_coeffs.span_cview());

    host_t<DFieldY> values(gridy);

    ddc::for_each(gridy, [&](ddc::DiscreteElement<IDimY> const idx) {
        values(idx) = sin(ddc::coordinate(idx));
    });
    double integral = integrate(values);
    return std::abs(2 - integral);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2)...};
}

TEST(SplineQuadratureTest, UniformConverge)
{
    constexpr int NTESTS(6);

    std::array<double, NTESTS> error = compute_errors(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        double order_error = abs(4 - order);
        EXPECT_LE(order_error, 5e-2);
    }
}

} // namespace
