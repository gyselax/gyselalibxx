// SPDX-License-Identifier: MIT


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "quadrature.hpp"
#include "simpson_quadrature.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

struct DimXPeriod
{
    static bool constexpr PERIODIC = true;
};

struct GridXPeriod : NonUniformGridBase<DimXPeriod>
{
};
using IDomXPeriod = IdxRange<GridXPeriod>;
using CoordXPeriod = Coord<DimXPeriod>;

TEST(TrapezoidUniformPeriodicQuadrature1D, ExactForConstantFunc)
{
    CoordXPeriod const x_min(0.0);
    CoordXPeriod const x_max(M_PI);
    IdxStep<GridXPeriod> const x_size(100);

    // Creating mesh & supports
    std::vector<CoordXPeriod> point_sampling;
    double dx = (x_max - x_min) / x_size;

    // Create & intialize Uniform mesh
    for (int k = 0; k < x_size; k++) {
        point_sampling.push_back(x_min + k * dx);
    }

    ddc::init_discrete_space<GridXPeriod>(point_sampling);
    Idx<GridXPeriod> lbound(0);
    IdxStep<GridXPeriod> npoints(x_size);
    IdxRange<GridXPeriod> gridx(lbound, npoints);

    host_t<FieldMem<double, IDomXPeriod>> const quadrature_coeffs_host(
            trapezoid_quadrature_coefficients(gridx));
    auto quadrature_coeffs = ddc::create_mirror_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(quadrature_coeffs_host));
    Quadrature const integrate(quadrature_coeffs.span_cview());

    FieldMem<double, IDomXPeriod> values(gridx);
    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), get_field(values), 1.0);
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values.span_cview());
    double expected_val = x_max - x_min;
    EXPECT_LE(abs(integral - expected_val), 1e-9);
}

TEST(SimpsonUniformPeriodicQuadrature1D, ExactForConstantFunc)
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

    host_t<FieldMem<double, IDomXPeriod>> const quadrature_coeffs_host(
            simpson_quadrature_coefficients_1d(gridx));
    auto quadrature_coeffs = ddc::create_mirror_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(quadrature_coeffs_host));
    Quadrature const integrate(quadrature_coeffs.span_cview());
    FieldMem<double, IDomXPeriod> values(gridx);

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), get_field(values), 1.0);
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values.span_cview());
    double expected_val = x_max - x_min;

    EXPECT_LE(abs(integral - expected_val), 1e-9);
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


enum Method { TRAPEZ, SIMPSON };

template <std::size_t N>
double compute_error(int n_elems, Method meth)
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

    host_t<DFieldMemY> quadrature_coeffs_host;
    switch (meth) {
    case Method::TRAPEZ:
        quadrature_coeffs_host = trapezoid_quadrature_coefficients(gridy);
    case Method::SIMPSON:
        quadrature_coeffs_host = simpson_quadrature_coefficients_1d(gridy);
    }

    auto quadrature_coeffs = ddc::create_mirror_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(quadrature_coeffs_host));

    Quadrature const integrate(quadrature_coeffs.span_cview());

    DFieldMemY values_alloc(gridy);
    DFieldY values = get_field(values_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridy,
            KOKKOS_LAMBDA(Idx<GridY> const idx) {
                values(idx) = Kokkos::sin(ddc::coordinate(idx));
            });
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values.span_cview());
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

TEST(TrapezoidUniformPeriodicQuadrature1D, Convergence)
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


TEST(SimpsonUniformPeriodicQuadrature1D, Convergence)
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
