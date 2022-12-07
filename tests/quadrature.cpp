// SPDX-License-Identifier: MIT

#include <sll/bsplines_non_uniform.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

TEST(QuadratureTest, ExactForConstantFunc)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IVectX const x_size(10);

    // Creating mesh & supports
    init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    init_discrete_space<IDimX>(InterpPointsX::get_sampling());
    DiscreteDomain<IDimX> interpolation_domain_x(InterpPointsX::get_domain());

    SplineXBuilder const builder_x(interpolation_domain_x);

    IDomainX const gridx = builder_x.interpolation_domain();

    Quadrature<IDimX> const integrate(trapezoid_quadrature_coefficients(gridx));

    DFieldX values(gridx);

    for_each(gridx, [&](DiscreteElement<IDimX> const idx) { values(idx) = 1.0; });
    double integral = integrate(values);
    double expected_val = x_max - x_min;
    EXPECT_LE(abs(integral - expected_val), 1e-9);
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
    using IDomainY = DiscreteDomain<IDimY>;
    using DFieldY = Chunk<double, IDomainY>;

    Coordinate<Y<N>> const x_min(0.0);
    Coordinate<Y<N>> const x_max(M_PI);

    init_discrete_space<BSplinesY>(x_min, x_max, n_elems);

    init_discrete_space<IDimY>(GrevillePointsY::get_sampling());
    IDomainY const gridx(GrevillePointsY::get_domain());

    SplineYBuilder const builder_x(gridx);

    Quadrature<IDimY> const integrate(trapezoid_quadrature_coefficients(gridx));

    DFieldY values(gridx);

    for_each(gridx, [&](DiscreteElement<IDimY> const idx) { values(idx) = sin(coordinate(idx)); });
    double integral = integrate(values);
    return std::abs(2 - integral);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2)...};
}

TEST(QuadratureTest, UniformConverge)
{
    constexpr int NTESTS(10);

    std::array<double, NTESTS> error = compute_errors(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        double order_error = abs(2 - order);
        EXPECT_LE(order_error, 1e-2);
    }
}
