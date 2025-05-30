// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "quadrature.hpp"
#include "simpson_quadrature.hpp"
#include "trapezoid_quadrature.hpp"


namespace {

enum Method { TRAPEZ, SIMPSON };

struct X
{
    static constexpr bool PERIODIC = false;
};

struct Y
{
    static constexpr bool PERIODIC = false;
};

using CoordX = Coord<X>;
using CoordY = Coord<Y>;

struct GridX : UniformGridBase<X>
{
};

struct GridY : UniformGridBase<Y>
{
};

using IdxXY = Idx<GridX, GridY>;

using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;

using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

using DFieldMemXY = DFieldMem<IdxRangeXY>;

double constant_func_check_2d(Method quad_method)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IdxStepX const x_ncells(10);

    CoordY const y_min(0.0);
    CoordY const y_max(20.0);
    IdxStepY const y_ncells(10);

    // Creating mesh & supports
    IdxRangeX const gridx
            = ddc::init_discrete_space<GridX>(GridX::init(x_min, x_max, x_ncells + 1));
    IdxRangeY const gridy
            = ddc::init_discrete_space<GridY>(GridY::init(y_min, y_max, y_ncells + 1));

    IdxRangeXY const gridxy(gridx, gridy);

    DFieldMemXY quadrature_coeffs;
    switch (quad_method) {
    case Method::TRAPEZ: {
        quadrature_coeffs
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);
        break;
    }
    case Method::SIMPSON: {
        quadrature_coeffs = simpson_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);
        break;
    }
    }

    Quadrature<IdxRangeXY> const integrate(get_const_field(quadrature_coeffs));
    DFieldMemXY values(gridxy);

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), values, 1.0);
    double integral = integrate(Kokkos::DefaultExecutionSpace(), get_const_field(values));
    double expected_val = (x_max - x_min) * (y_max - y_min);

    return abs(integral - expected_val);
}

void integrated_function_operator(Method quad_method)
{
    CoordX x_min(0.0);
    CoordX x_max(3.0);
    IdxStepX x_ncells(4);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IdxStepY y_ncells(16);

    IdxRangeX gridx
            = ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_ncells + 1));
    IdxRangeY gridy
            = ddc::init_discrete_space<GridY>(GridY::init<GridY>(y_min, y_max, y_ncells + 1));
    IdxRangeXY gridxy(gridx, gridy);

    DFieldMemXY quadrature_coeffs;
    switch (quad_method) {
    case Method::TRAPEZ: {
        quadrature_coeffs
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);
        break;
    }
    case Method::SIMPSON: {
        quadrature_coeffs = simpson_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);
        break;
    }
    }
    Quadrature func_operator(get_const_field(quadrature_coeffs));

    double const integral = func_operator(
            Kokkos::DefaultExecutionSpace(),
            KOKKOS_LAMBDA(IdxXY ixy) {
                double y = ddc::coordinate(ddc::select<GridY>(ixy));
                double x = ddc::coordinate(ddc::select<GridX>(ixy));
                return x * y + 2;
            });
    EXPECT_DOUBLE_EQ(integral, 132.);
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
    struct GridX : UniformGridBase<X>
    {
    };
    struct GridY : UniformGridBase<Y>
    {
    };
};

template <std::size_t N>
double compute_error(int n_elems, Method quad_method)
{
    using DimX = typename ComputeErrorTraits<N>::X;
    using DimY = typename ComputeErrorTraits<N>::Y;
    using GridX = typename ComputeErrorTraits<N>::GridX;
    using GridY = typename ComputeErrorTraits<N>::GridY;
    using IdxStepX = IdxStep<GridX>;
    using IdxStepY = IdxStep<GridY>;
    using IdxRangeX = IdxRange<GridX>;
    using IdxRangeY = IdxRange<GridY>;
    using IdxRangeXY = IdxRange<GridX, GridY>;
    using DFieldMemXY = DFieldMem<IdxRangeXY>;
    using DFieldXY = DField<IdxRangeXY>;

    Coord<DimX> const x_min(0.0);
    Coord<DimX> const x_max(M_PI);
    IdxStepX x_ncells(n_elems);

    Coord<DimY> const y_min(0.0);
    Coord<DimY> const y_max(M_PI);
    IdxStepY y_ncells(n_elems);

    IdxRangeX const gridx = ddc::init_discrete_space<GridX>(
            GridX::template init<GridX>(x_min, x_max, x_ncells + 1));
    IdxRangeY const gridy = ddc::init_discrete_space<GridY>(
            GridY::template init<GridY>(y_min, y_max, y_ncells + 1));
    IdxRangeXY const gridxy(gridx, gridy);

    DFieldMemXY quadrature_coeffs;
    switch (quad_method) {
    case Method::TRAPEZ: {
        quadrature_coeffs
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);
        break;
    }
    case Method::SIMPSON: {
        quadrature_coeffs = simpson_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);
        break;
    }
    }
    Quadrature<IdxRangeXY> const integrate(get_const_field(quadrature_coeffs));

    DFieldMemXY values_alloc(gridxy);
    DFieldXY values = get_field(values_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(Idx<GridX, GridY> const idx) {
                double const y_cos = Kokkos::cos(ddc::get<DimY>(ddc::coordinate(idx)));
                values(idx) = sin(ddc::get<DimX>(ddc::coordinate(idx))) * y_cos * y_cos;
            });
    double integral = integrate(Kokkos::DefaultExecutionSpace(), values);
    return std::abs(integral - M_PI);
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_trapez_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2, Method::TRAPEZ)...};
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_simps_errors(std::index_sequence<Is...>, int n_elems)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(n_elems *= 2, Method::SIMPSON)...};
}

TEST(TrapezoidUniformNonPeriodicQuadrature, ExactForLinearBatchSecond2D)
{
    integrated_function_operator(Method::TRAPEZ);
}

TEST(SimpsonUniformNonPeriodicQuadrature, ExactForLinearBatchSecond2D)
{
    integrated_function_operator(Method::SIMPSON);
}

TEST(TrapezoidUniformNonPeriodicQuadrature2D, ExactForConstantFunc)
{
    EXPECT_LE(constant_func_check_2d(Method::TRAPEZ), 1e-9);
}

TEST(SimpsonUniformNonPeriodicQuadrature2D, ExactForConstantFunc)
{
    EXPECT_LE(constant_func_check_2d(Method::SIMPSON), 1e-9);
}

TEST(TrapezoidUniformNonPeriodicQuadrature2D, Convergence)
{
    constexpr int NTESTS(4);

    std::array<double, NTESTS> error
            = compute_trapez_errors(std::make_index_sequence<NTESTS>(), 50);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        EXPECT_NEAR(order, 2, 1e-1);
    }
}

TEST(SimpsonUniformNonPeriodicQuadrature2D, Convergence)
{
    constexpr int NTESTS(4);

    std::array<double, NTESTS> error = compute_simps_errors(std::make_index_sequence<NTESTS>(), 50);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        EXPECT_NEAR(order, 4, 1e-1);
    }
}
} // namespace
