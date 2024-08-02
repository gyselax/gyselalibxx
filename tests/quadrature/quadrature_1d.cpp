// SPDX-License-Identifier: MIT


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "quadrature.hpp"
#include "simpson_quadrature.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

enum Method { TRAPEZ, SIMPSON };

struct DimXPeriod
{
    static bool constexpr PERIODIC = true;
};

struct GridXPeriod : NonUniformGridBase<DimXPeriod>
{
};
using IDomXPeriod = IdxRange<GridXPeriod>;
using CoordXPeriod = Coord<DimXPeriod>;
using DFieldMemX = FieldMem<double, IDomXPeriod>;

double constant_func_check_1d(Method quad_method)
{
    CoordXPeriod const x_min(0.0);
    CoordXPeriod const x_max(M_PI);
    IdxStep<GridXPeriod> const x_size(100);

    // Creating mesh & support
    std::vector<CoordXPeriod> point_sampling;
    double dx = (x_max - x_min) / x_size;
    // Create & intialize mesh
    for (int k = 0; k < x_size; k++) {
        point_sampling.push_back(x_min + k * dx);
    }

    ddc::init_discrete_space<GridXPeriod>(point_sampling);
    Idx<GridXPeriod> lbound(0);
    IdxStep<GridXPeriod> npoints(x_size);
    IdxRange<GridXPeriod> gridx(lbound, npoints);

    DFieldMemX quadrature_coeffs_alloc;
    switch (quad_method) {
    case Method::TRAPEZ:
        quadrature_coeffs_alloc
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx);
        break;
    case Method::SIMPSON:
        quadrature_coeffs_alloc
                = simpson_quadrature_coefficients_1d<Kokkos::DefaultExecutionSpace>(gridx);
        break;
    }

    Quadrature const integrate(get_const_field(quadrature_coeffs_alloc));

    DFieldMemX values_alloc(gridx);
    device_t<Field<double, IDomXPeriod>> values = get_field(values_alloc);
    ddc::parallel_fill(values, 1.0);
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values);
    double expected_val = x_max - x_min;

    return abs(integral - expected_val);
}

template <std::size_t N>
struct ComputeErrorTraits
{
    struct Y
    {
        static bool constexpr PERIODIC = false;
    };
    struct GridY : NonUniformGridBase<Y>
    {
    };
};

template <std::size_t N>
double compute_error(int n_elems, Method quad_method)
{
    using DimY = typename ComputeErrorTraits<N>::Y;
    using GridY = typename ComputeErrorTraits<N>::GridY;
    using IdxRangeY = IdxRange<GridY>;
    using DFieldMemY = FieldMem<double, IdxRangeY>;
    using DFieldY = Field<double, IdxRangeY>;

    Coord<DimY> const y_min(0.0);
    Coord<DimY> const y_max(M_PI);
    std::vector<Coord<DimY>> point_sampling;
    double dy = (y_max - y_min) / n_elems;

    // Create & intialize Uniform mesh
    for (int k = 0; k < n_elems; k++) {
        point_sampling.push_back(y_min + k * dy);
    }

    ddc::init_discrete_space<GridY>(point_sampling);
    Idx<GridY> lbound(0);
    IdxStep<GridY> npoints(n_elems);
    IdxRange<GridY> gridy(lbound, npoints);

    DFieldMemX quadrature_coeffs_alloc;
    switch (quad_method) {
    case Method::TRAPEZ:
        quadrature_coeffs_alloc
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridy);
        break;
    case Method::SIMPSON:
        quadrature_coeffs_alloc
                = simpson_quadrature_coefficients_1d<Kokkos::DefaultExecutionSpace>(gridy);
        break;
    }

    Quadrature const integrate(get_const_field(quadrature_coeffs_alloc));
    DFieldMemY values_alloc(gridy);
    DFieldY values = get_field(values_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridy,
            KOKKOS_LAMBDA(Idx<GridY> const idx) { values(idx) = sin(ddc::coordinate(idx)); });

    double integral = integrate(Kokkos::DefaultExecutionSpace(), values);
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

TEST(TrapezoidUniformPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(constant_func_check_1d(Method::TRAPEZ), 1e-9);
}

TEST(SimpsonUniformPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(constant_func_check_1d(Method::SIMPSON), 1e-9);
}

TEST(TrapezoidUniformNonPeriodicQuadrature1D, Convergence)
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

TEST(SimpsonUniformNonPeriodicQuadrature1D, Convergence)
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
