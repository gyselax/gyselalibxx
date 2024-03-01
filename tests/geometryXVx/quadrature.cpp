// SPDX-License-Identifier: MIT


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "quadrature.hpp"
#include "simpson_quadrature.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

struct RDimXPeriod
{
    static bool constexpr PERIODIC = true;
};

using IDimXPeriod = ddc::NonUniformPointSampling<RDimXPeriod>;
using IDomXPeriod = ddc::DiscreteDomain<IDimXPeriod>;
using CoordXPeriod = ddc::Coordinate<RDimXPeriod>;

TEST(QuadratureTest, TrapezExactForConstantFunc)
{
    CoordXPeriod const x_min(0.0);
    CoordXPeriod const x_max(M_PI);
    ddc::DiscreteVector<IDimXPeriod> const x_size(100);

    // Creating mesh & supports
    std::vector<CoordXPeriod> point_sampling;
    double dx = (x_max - x_min) / x_size;

    // Create & intialize Uniform mesh
    for (int k = 0; k < x_size; k++) {
        point_sampling.push_back(x_min + k * dx);
    }

    ddc::init_discrete_space<IDimXPeriod>(point_sampling);
    ddc::DiscreteElement<IDimXPeriod> lbound(0);
    ddc::DiscreteVector<IDimXPeriod> npoints(x_size);
    ddc::DiscreteDomain<IDimXPeriod> gridx(lbound, npoints);

    ddc::Chunk<double, IDomXPeriod> const quadrature_coeffs
            = trapezoid_quadrature_coefficients(gridx);
    Quadrature<IDimXPeriod> const integrate(quadrature_coeffs);

    ddc::Chunk<double, IDomXPeriod> values(gridx);

    ddc::for_each(gridx, [&](ddc::DiscreteElement<IDimXPeriod> const idx) { values(idx) = 1.0; });
    double integral = integrate(values);
    double expected_val = x_max - x_min;
    EXPECT_LE(abs(integral - expected_val), 1e-9);
}

TEST(QuadratureTest, SimpsonExactForConstantFunc)
{
    CoordXPeriod const x_min(0.0);
    CoordXPeriod const x_max(M_PI);
    ddc::DiscreteVector<IDimXPeriod> const x_size(100);

    // Creating mesh & support
    std::vector<CoordXPeriod> point_sampling;
    double dx = (x_max - x_min) / x_size;
    // Create & intialize mesh
    for (int k = 0; k < x_size; k++) {
        point_sampling.push_back(x_min + k * dx);
    }

    ddc::init_discrete_space<IDimXPeriod>(point_sampling);
    ddc::DiscreteElement<IDimXPeriod> lbound(0);
    ddc::DiscreteVector<IDimXPeriod> npoints(x_size);
    ddc::DiscreteDomain<IDimXPeriod> gridx(lbound, npoints);

    ddc::Chunk<double, IDomXPeriod> const quadrature_coeffs
            = simpson_quadrature_coefficients_1d(gridx);
    Quadrature<IDimXPeriod> const integrate(quadrature_coeffs);
    ddc::Chunk<double, IDomXPeriod> values(gridx);

    ddc::for_each(gridx, [&](ddc::DiscreteElement<IDimXPeriod> const idx) { values(idx) = 1.0; });
    double integral = integrate(values);
    double expected_val = x_max - x_min;

    EXPECT_LE(abs(integral - expected_val), 1e-9);
}

template <std::size_t N>
struct Y
{
    static bool constexpr PERIODIC = false;
};


enum Method { TRAPEZ, SIMPSON };

template <std::size_t N>
double compute_error(int n_elems, Method meth)
{
    using DimY = Y<N>;
    using IDimY = ddc::NonUniformPointSampling<DimY>;
    using IDomainY = ddc::DiscreteDomain<IDimY>;
    using DFieldY = device_t<ddc::Chunk<double, IDomainY>>;

    ddc::Coordinate<DimY> const y_min(0.0);
    ddc::Coordinate<DimY> const y_max(M_PI);
    std::vector<ddc::Coordinate<DimY>> point_sampling;
    double dy = (y_max - y_min) / n_elems;

    // Create & intialize Uniform mesh
    for (int k = 0; k < n_elems; k++) {
        point_sampling.push_back(y_min + k * dy);
    }

    ddc::init_discrete_space<IDimY>(point_sampling);
    ddc::DiscreteElement<IDimY> lbound(0);
    ddc::DiscreteVector<IDimY> npoints(n_elems);
    ddc::DiscreteDomain<IDimY> gridy(lbound, npoints);

    host_t<DFieldY> quadrature_coeffs;
    switch (meth) {
    case Method::TRAPEZ:
        quadrature_coeffs = trapezoid_quadrature_coefficients(gridy);
    case Method::SIMPSON:
        quadrature_coeffs = simpson_quadrature_coefficients_1d(gridy);
    }

    Quadrature<IDimY> const integrate(quadrature_coeffs);
    host_t<DFieldY> values(gridy);

    ddc::for_each(gridy, [&](ddc::DiscreteElement<IDimY> const idx) {
        values(idx) = sin(ddc::coordinate(idx));
    });
    double integral = integrate(values);
    return std::abs(2 - integral);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors_trpz(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2, Method::TRAPEZ)...};
}
template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors_simpson(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2, Method::SIMPSON)...};
}

TEST(QuadratureTest, TrapezUniformConverge)
{
    constexpr int NTESTS(10);

    std::array<double, NTESTS> error = compute_errors_trpz(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        double order_error = abs(2 - order);
        EXPECT_LE(order_error, 1e-2);
    }
}


TEST(QuadratureTest, SimpsonUniformConverge)
{
    constexpr int NTESTS(10);

    std::array<double, NTESTS> error
            = compute_errors_simpson(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        double order_error = abs(2 - order);
        EXPECT_LE(order_error, 1e-2);
    }
}

} // namespace
