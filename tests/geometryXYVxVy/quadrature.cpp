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

    CoordY const y_min(0.0);
    CoordY const y_max(20.0);
    IVectY const y_size(10);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesY>(y_min, y_max, y_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimY>(SplineInterpPointsY::get_sampling());

    IDomainX const gridx(SplineInterpPointsX::get_domain());
    IDomainY const gridy(SplineInterpPointsY::get_domain());

    IDomainXY const gridxy(gridx, gridy);

    Quadrature<IDimX, IDimY> const integrate(trapezoid_quadrature_coefficients(gridxy));

    DFieldXY values(gridxy);

    ddc::for_each(gridxy, [&](ddc::DiscreteElement<IDimX, IDimY> const idx) { values(idx) = 1.0; });
    double integral = integrate(values);
    double expected_val = (x_max - x_min) * (y_max - y_min);
    EXPECT_LE(abs(integral - expected_val), 1e-9);
}

template <std::size_t N>
struct X
{
    static bool constexpr PERIODIC = false;
};

template <std::size_t N>
struct Y
{
    static bool constexpr PERIODIC = false;
};

template <std::size_t N>
double compute_error(int n_elems)
{
    using DimX = X<N>;
    using DimY = Y<N>;
    using BSplinesX = UniformBSplines<DimX, 3>;
    using BSplinesY = UniformBSplines<DimY, 3>;
    using GrevillePointsX
            = GrevilleInterpolationPoints<BSplinesX, BoundCond::GREVILLE, BoundCond::GREVILLE>;
    using GrevillePointsY
            = GrevilleInterpolationPoints<BSplinesY, BoundCond::GREVILLE, BoundCond::GREVILLE>;
    using IDimX = typename GrevillePointsX::interpolation_mesh_type;
    using IDimY = typename GrevillePointsY::interpolation_mesh_type;
    using IDomainX = ddc::DiscreteDomain<IDimX>;
    using IDomainY = ddc::DiscreteDomain<IDimY>;
    using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
    using DFieldXY = ddc::Chunk<double, IDomainXY>;

    ddc::Coordinate<DimX> const x_min(0.0);
    ddc::Coordinate<DimX> const x_max(M_PI);

    ddc::Coordinate<DimY> const y_min(0.0);
    ddc::Coordinate<DimY> const y_max(M_PI);

    ddc::init_discrete_space<BSplinesX>(x_min, x_max, n_elems);
    ddc::init_discrete_space<BSplinesY>(y_min, y_max, n_elems);

    ddc::init_discrete_space<IDimX>(GrevillePointsX::get_sampling());
    ddc::init_discrete_space<IDimY>(GrevillePointsY::get_sampling());
    IDomainX const gridx(GrevillePointsX::get_domain());
    IDomainY const gridy(GrevillePointsY::get_domain());
    IDomainXY const gridxy(gridx, gridy);

    Quadrature<IDimX, IDimY> const integrate(trapezoid_quadrature_coefficients(gridxy));

    DFieldXY values(gridxy);

    ddc::for_each(gridxy, [&](ddc::DiscreteElement<IDimX, IDimY> const idx) {
        double const y_cos = cos(ddc::get<DimY>(ddc::coordinate(idx)));
        values(idx) = sin(ddc::get<DimX>(ddc::coordinate(idx))) * y_cos * y_cos;
    });
    double integral = integrate(values);
    return std::abs(integral - M_PI);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2)...};
}

TEST(QuadratureTest, UniformConverge)
{
    constexpr int NTESTS(4);

    std::array<double, NTESTS> error = compute_errors(std::make_index_sequence<NTESTS>(), 50);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        double order_error = abs(2 - order);
        EXPECT_LE(order_error, 1e-1);
    }
}
