// SPDX-License-Identifier: MIT


#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "quadrature.hpp"
#include "simpson_quadrature.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

enum Method { TRAPEZ, SIMPSON };

struct RDimXPeriod
{
    static bool constexpr PERIODIC = true;
};

struct IDimXPeriod : ddc::NonUniformPointSampling<RDimXPeriod>
{
};
using IDomXPeriod = ddc::DiscreteDomain<IDimXPeriod>;
using CoordXPeriod = ddc::Coordinate<RDimXPeriod>;
using DFieldX = device_t<ddc::Chunk<double, IDomXPeriod>>;

double constant_func_check_1d(Method quad_method)
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

    DFieldX quadrature_coeffs_alloc(gridx);
    if (quad_method == Method::TRAPEZ) {
        quadrature_coeffs_alloc
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx);
    }
    if (quad_method == Method::SIMPSON) {
        quadrature_coeffs_alloc
                = simpson_quadrature_coefficients_1d<Kokkos::DefaultExecutionSpace>(gridx);
    }

    Quadrature const integrate(quadrature_coeffs_alloc.span_cview());

    DFieldX values_alloc(gridx);
    ddc::ChunkSpan values = values_alloc.span_view();
    Kokkos::deep_copy(values.allocation_kokkos_view(), 1.0);
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
    struct IDimY : ddc::NonUniformPointSampling<Y>
    {
    };
};

template <std::size_t N>
double compute_error(int n_elems, Method quad_method)
{
    using DimY = typename ComputeErrorTraits<N>::Y;
    using IDimY = typename ComputeErrorTraits<N>::IDimY;
    using IDomainY = ddc::DiscreteDomain<IDimY>;
    using DFieldY = device_t<ddc::Chunk<double, IDomainY>>;
    using DSpanY = device_t<ddc::ChunkSpan<double, IDomainY>>;

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

    DFieldY quadrature_coeffs_alloc(gridy);
    if (quad_method == Method::TRAPEZ) {
        quadrature_coeffs_alloc
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace, IDimY>(gridy);
    }
    if (quad_method == Method::SIMPSON) {
        quadrature_coeffs_alloc
                = simpson_quadrature_coefficients_1d<Kokkos::DefaultExecutionSpace, IDimY>(gridy);
    }

    Quadrature const integrate(quadrature_coeffs_alloc.span_cview());
    DFieldY values_alloc(gridy);
    ddc::ChunkSpan values = values_alloc.span_view();

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridy,
            KOKKOS_LAMBDA(ddc::DiscreteElement<IDimY> const idx) {
                values(idx) = sin(ddc::coordinate(idx));
            });

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
