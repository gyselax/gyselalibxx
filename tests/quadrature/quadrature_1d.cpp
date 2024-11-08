// SPDX-License-Identifier: MIT
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "mesh_builder.hpp"
#include "quadrature.hpp"
#include "simpson_quadrature.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

enum Method { TRAPEZ, SIMPSON };

template <bool Periodic>
struct ConstantFuncCheck
{
    struct XPeriod
    {
        static bool constexpr PERIODIC = Periodic;
    };

    struct GridXPeriod : NonUniformGridBase<XPeriod>
    {
    };
    using IdxRangeXPeriod = IdxRange<GridXPeriod>;
    using CoordXPeriod = Coord<XPeriod>;
    using DFieldMemX = DFieldMem<IdxRangeXPeriod>;

    static double check_1d(Method quad_method)
    {
        CoordXPeriod const x_min(0.0);
        CoordXPeriod const x_max(1.0);
        IdxStep<GridXPeriod> const x_size(100);

        // Creating mesh & support
        std::vector<CoordXPeriod> point_sampling
                = build_random_non_uniform_break_points(x_min, x_max, x_size);
        IdxStep<GridXPeriod> npoints(x_size + 1);
        ddc::init_discrete_space<GridXPeriod>(point_sampling);
        Idx<GridXPeriod> lbound(0);
        IdxRange<GridXPeriod> gridx(lbound, npoints - int(Periodic));

        DFieldMemX quadrature_coeffs_alloc;
        switch (quad_method) {
        case Method::TRAPEZ: {
            quadrature_coeffs_alloc
                    = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx);
            break;
        }
        case Method::SIMPSON: {
            quadrature_coeffs_alloc
                    = simpson_quadrature_coefficients_1d<Kokkos::DefaultExecutionSpace>(gridx);
            break;
        }
        }

        Quadrature const integrate(get_const_field(quadrature_coeffs_alloc));

        DFieldMemX values_alloc(gridx);
        DField<IdxRangeXPeriod> values = get_field(values_alloc);
        ddc::parallel_fill(values, 1.0);
        double integral = integrate(Kokkos::DefaultExecutionSpace(), values);
        double expected_val = x_max - x_min;

        return abs(integral - expected_val);
    }

    static double check_simpson_trapezoid_1d(Extremity trapezoid_extremity)
    {
        CoordXPeriod const x_min(0.0);
        CoordXPeriod const x_max(1.0);
        IdxStep<GridXPeriod> const x_size(5);

        // Creating mesh & support
        std::vector<CoordXPeriod> point_sampling
                = build_random_non_uniform_break_points(x_min, x_max, x_size);
        ddc::init_discrete_space<GridXPeriod>(point_sampling);
        Idx<GridXPeriod> lbound(0);
        IdxStep<GridXPeriod> npoints(x_size + 1);
        IdxRange<GridXPeriod> gridx(lbound, npoints);

        DFieldMemX quadrature_coeffs_alloc = simpson_trapezoid_quadrature_coefficients_1d<
                Kokkos::DefaultExecutionSpace>(gridx, trapezoid_extremity);

        Quadrature const integrate(get_const_field(quadrature_coeffs_alloc));

        DFieldMemX values(gridx);
        ddc::parallel_fill(get_field(values), 1.0);
        double integral = integrate(Kokkos::DefaultExecutionSpace(), get_const_field(values));
        double expected_val = x_max - x_min;
        return abs(integral - expected_val);
    }
};

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
double compute_error(int n_cells, Method quad_method)
{
    using DimY = typename ComputeErrorTraits<N>::Y;
    using GridY = typename ComputeErrorTraits<N>::GridY;
    using IdxRangeY = IdxRange<GridY>;
    using DFieldMemY = DFieldMem<IdxRangeY>;
    using DFieldY = DField<IdxRangeY>;

    Coord<DimY> const y_min(0.0);
    Coord<DimY> const y_max(M_PI);
    IdxStep<GridY> npoints(n_cells + 1);
    IdxStep<GridY> ncells(n_cells);

    std::vector<Coord<DimY>> point_sampling = build_uniform_break_points(y_min, y_max, ncells);
    ddc::init_discrete_space<GridY>(point_sampling);

    Idx<GridY> lbound(0);
    IdxRange<GridY> gridy(lbound, npoints);

    DFieldMemY quadrature_coeffs_alloc;
    switch (quad_method) {
    case Method::TRAPEZ: {
        quadrature_coeffs_alloc
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridy);
        break;
    }
    case Method::SIMPSON: {
        quadrature_coeffs_alloc
                = simpson_quadrature_coefficients_1d<Kokkos::DefaultExecutionSpace>(gridy);
        break;
    }
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
std::array<double, sizeof...(Is)> compute_errors_trpz(std::index_sequence<Is...>, int ncells)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(ncells *= 2, Method::TRAPEZ)...};
}

template <std::size_t... Is>
std::array<double, sizeof...(Is)> compute_errors_simpson(std::index_sequence<Is...>, int ncells)
{
    return std::array<double, sizeof...(Is)> {compute_error<Is>(ncells *= 2, Method::SIMPSON)...};
}

TEST(TrapezoidNonUniformPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(ConstantFuncCheck<true>::check_1d(Method::TRAPEZ), 1e-11);
}

TEST(SimpsonNonUniformPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(ConstantFuncCheck<true>::check_1d(Method::SIMPSON), 1e-11);
}

TEST(TrapezoidNonUniformNonPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(ConstantFuncCheck<false>::check_1d(Method::TRAPEZ), 1e-11);
}

TEST(SimpsonNonUniformNonPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(ConstantFuncCheck<false>::check_1d(Method::SIMPSON), 1e-11);
}

TEST(SimpsonTrapezoidFrontNonUniformNonPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(ConstantFuncCheck<false>::check_simpson_trapezoid_1d(Extremity::FRONT), 1e-11);
}

TEST(SimpsonTrapezoidBackNonUniformNonPeriodicQuadrature1D, ExactForConstantFunc)
{
    EXPECT_LE(ConstantFuncCheck<false>::check_simpson_trapezoid_1d(Extremity::BACK), 1e-11);
}

TEST(TrapezoidUniformNonPeriodicQuadrature1D, Convergence)
{
    constexpr int NTESTS(10);

    std::array<double, NTESTS> error = compute_errors_trpz(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        EXPECT_NEAR(order, 2, 1e-2);
    }
}

TEST(SimpsonUniformNonPeriodicQuadrature1D, Convergence)
{
    constexpr int NTESTS(5);

    std::array<double, NTESTS> error
            = compute_errors_simpson(std::make_index_sequence<NTESTS>(), 10);

    for (int i(1); i < NTESTS; ++i) {
        EXPECT_LE(error[i], error[i - 1]);
        double order = log(error[i - 1] / error[i]) / log(2.0);
        EXPECT_NEAR(order, 4, 1e-2);
    }
}

} // namespace
