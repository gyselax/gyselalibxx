// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_helper.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

struct X
{
    static constexpr bool PERIODIC = false;
};

struct Y
{
    static constexpr bool PERIODIC = false;
};

using CoordX = ddc::Coordinate<X>;
using CoordY = ddc::Coordinate<Y>;

struct IDimX : ddc::UniformPointSampling<X>
{
};

struct IDimY : ddc::UniformPointSampling<Y>
{
};

using IdxXY = ddc::DiscreteElement<IDimX, IDimY>;

using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;

using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;

using DFieldXY = device_t<ddc::Chunk<double, IDomainXY>>;

TEST(TrapezoidUniformNonPeriodicQuadrature2D, ExactForConstantFunc)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IVectX const x_size(10);

    CoordY const y_min(0.0);
    CoordY const y_max(20.0);
    IVectY const y_size(10);

    // Creating mesh & supports
    IDomainX const gridx = ddc::init_discrete_space<IDimX>(IDimX::init(x_min, x_max, x_size));
    IDomainY const gridy = ddc::init_discrete_space<IDimY>(IDimY::init(y_min, y_max, y_size));

    IDomainXY const gridxy(gridx, gridy);

    host_t<DFieldXY> const quadrature_coeffs_host = trapezoid_quadrature_coefficients(gridxy);
    auto quadrature_coeffs = ddc::create_mirror_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    Quadrature<IDomainXY> const integrate(quadrature_coeffs);

    DFieldXY values(gridxy);

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), values, 1.0);
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values.span_cview());
    double expected_val = (x_max - x_min) * (y_max - y_min);
    EXPECT_LE(abs(integral - expected_val), 1e-9);
}

void integrated_function_operator()
{
    CoordX x_min(0.0);
    CoordX x_max(3.0);
    IVectX x_ncells(4);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IVectY y_ncells(16);

    IDomainX gridx = ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_ncells));

    IDomainY gridy = ddc::init_discrete_space<IDimY>(IDimY::init<IDimY>(y_min, y_max, y_ncells));
    IDomainXY gridxy(gridx, gridy);

    host_t<DFieldXY> quad_coeffs_host_second = trapezoid_quadrature_coefficients(gridxy);
    auto quad_coeffs_second = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quad_coeffs_host_second.span_view());
    Quadrature func_operator(quad_coeffs_second.span_cview());

    double const integral = func_operator(
            Kokkos::DefaultExecutionSpace(),
            KOKKOS_LAMBDA(IdxXY ixy) {
                double y = ddc::coordinate(ddc::select<IDimY>(ixy));
                double x = ddc::coordinate(ddc::select<IDimX>(ixy));
                return x * y + 2;
            });
    EXPECT_DOUBLE_EQ(integral, 132.);
}

TEST(TrapezoidUniformNonPeriodicQuadrature, ExactForLinearBatchSecond2D)
{
    integrated_function_operator();
}

template <std::size_t N>
struct ComputeErrorTraits
{
    struct X
    {
        static bool constexpr PERIODIC = false;
    };
    struct Y
    {
        static bool constexpr PERIODIC = false;
    };
    struct IDimX : ddc::UniformPointSampling<X>
    {
    };
    struct IDimY : ddc::UniformPointSampling<Y>
    {
    };
};

template <std::size_t N>
double compute_error(int n_elems)
{
    using DimX = typename ComputeErrorTraits<N>::X;
    using DimY = typename ComputeErrorTraits<N>::Y;
    using IDimX = typename ComputeErrorTraits<N>::IDimX;
    using IDimY = typename ComputeErrorTraits<N>::IDimY;
    using IVectX = ddc::DiscreteVector<IDimX>;
    using IVectY = ddc::DiscreteVector<IDimY>;
    using IDomainX = ddc::DiscreteDomain<IDimX>;
    using IDomainY = ddc::DiscreteDomain<IDimY>;
    using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
    using DFieldXY = device_t<ddc::Chunk<double, IDomainXY>>;
    using DSpanXY = device_t<ddc::ChunkSpan<double, IDomainXY>>;

    ddc::Coordinate<DimX> const x_min(0.0);
    ddc::Coordinate<DimX> const x_max(M_PI);
    IVectX x_size(n_elems);

    ddc::Coordinate<DimY> const y_min(0.0);
    ddc::Coordinate<DimY> const y_max(M_PI);
    IVectY y_size(n_elems);

    IDomainX const gridx
            = ddc::init_discrete_space<IDimX>(IDimX::template init<IDimX>(x_min, x_max, x_size));
    IDomainY const gridy
            = ddc::init_discrete_space<IDimY>(IDimY::template init<IDimY>(y_min, y_max, y_size));
    IDomainXY const gridxy(gridx, gridy);

    host_t<DFieldXY> const quadrature_coeffs_host = trapezoid_quadrature_coefficients(gridxy);
    auto quadrature_coeffs = ddc::create_mirror_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    Quadrature<IDomainXY> const integrate(quadrature_coeffs);

    DFieldXY values_alloc(gridxy);
    DSpanXY values = values_alloc.span_view();

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(ddc::DiscreteElement<IDimX, IDimY> const idx) {
                double const y_cos = Kokkos::cos(ddc::get<DimY>(ddc::coordinate(idx)));
                values(idx) = sin(ddc::get<DimX>(ddc::coordinate(idx))) * y_cos * y_cos;
            });
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values);
    return std::abs(integral - M_PI);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2)...};
}

TEST(TrapezoidUniformNonPeriodicQuadrature2D, Convergence)
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
} // namespace
