// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_helper.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

namespace {

struct Batch1
{
};

struct Batch2
{
};

struct X
{
    static bool constexpr PERIODIC = false;
};

struct Y
{
    static bool constexpr PERIODIC = false;
};

using CoordBatch1 = const ddc::Coordinate<Batch1>;
using CoordBatch2 = const ddc::Coordinate<Batch2>;
using CoordX = const ddc::Coordinate<X>;
using CoordY = const ddc::Coordinate<Y>;

struct IDimBatch1 : ddc::UniformPointSampling<Batch1>
{
};
struct IDimBatch2 : ddc::UniformPointSampling<Batch2>
{
};

struct IDimX : ddc::UniformPointSampling<X>
{
};
struct IDimY : ddc::UniformPointSampling<Y>
{
};

using IdxBatch1 = ddc::DiscreteElement<IDimBatch1>;
using IdxBatch2 = ddc::DiscreteElement<IDimBatch2>;
using IdxX = ddc::DiscreteElement<IDimX>;
using IdxY = ddc::DiscreteElement<IDimY>;
using IdxXY = ddc::DiscreteElement<IDimX, IDimY>;
using IdxB1X = ddc::DiscreteElement<IDimBatch1, IDimX>;
using IdxB1B2 = ddc::DiscreteElement<IDimBatch1, IDimBatch2>;
using IdxB1XY = ddc::DiscreteElement<IDimBatch1, IDimX, IDimY>;
using IdxB1B2X = ddc::DiscreteElement<IDimBatch1, IDimBatch2, IDimX>;
using IdxB1B2XY = ddc::DiscreteElement<IDimBatch1, IDimBatch2, IDimX, IDimY>;
using IdxXB1YB2 = ddc::DiscreteElement<IDimX, IDimBatch1, IDimY, IDimBatch2>;

using IVectBatch1 = const ddc::DiscreteVector<IDimBatch1>;
using IVectBatch2 = const ddc::DiscreteVector<IDimBatch2>;
using IVectX = const ddc::DiscreteVector<IDimX>;
using IVectY = const ddc::DiscreteVector<IDimY>;

using IDomainBatch1 = ddc::DiscreteDomain<IDimBatch1>;
using IDomainBatch2 = ddc::DiscreteDomain<IDimBatch2>;
using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainB1X = ddc::DiscreteDomain<IDimBatch1, IDimX>;
using IDomainB1XY = ddc::DiscreteDomain<IDimBatch1, IDimX, IDimY>;
using IDomainB1B2X = ddc::DiscreteDomain<IDimBatch1, IDimBatch2, IDimX>;
using IDomainB1B2XY = ddc::DiscreteDomain<IDimBatch1, IDimBatch2, IDimX, IDimY>;
using IDomainXB1YB2 = ddc::DiscreteDomain<IDimX, IDimBatch1, IDimY, IDimBatch2>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainB1B2 = ddc::DiscreteDomain<IDimBatch1, IDimBatch2>;

using DFieldBatch1 = device_t<ddc::Chunk<double, IDomainBatch1>>;
using DFieldX = device_t<ddc::Chunk<double, IDomainX>>;
using DFieldY = device_t<ddc::Chunk<double, IDomainY>>;
using DFieldXY = device_t<ddc::Chunk<double, IDomainXY>>;
using DFieldB1B2 = device_t<ddc::Chunk<double, IDomainB1B2>>;

void batched_operator_1d()
{
    CoordBatch1 b_min(0.0);
    CoordBatch1 b_max(3.0);
    IVectBatch1 b_ncells(4);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IVectX x_ncells(16);

    IDomainBatch1 gridb = ddc::init_discrete_space<IDimBatch1>(
            IDimBatch1::init<IDimBatch1>(b_min, b_max, b_ncells));
    IDomainX gridx = ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_ncells));

    DFieldX quad_coeffs = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx);

    Quadrature<IDomainX, IDomainB1X> quad_batched_operator(quad_coeffs.span_view());

    DFieldBatch1 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            results.span_view(),
            KOKKOS_LAMBDA(IdxB1X ibx) {
                double b = ddc::coordinate(ddc::select<IDimBatch1>(ibx));
                double x = ddc::coordinate(ddc::select<IDimX>(ibx));
                return b * x + 2;
            });

    auto results_host = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), results.span_view());

    ddc::for_each(gridb, [&](IdxBatch1 ib) {
        double b = ddc::coordinate(ddc::select<IDimBatch1>(ib));
        double x = x_max;
        double const ubound = 0.5 * b * x * x + 2 * x;
        x = x_min;
        double const lbound = 0.5 * b * x * x + 2 * x;
        EXPECT_DOUBLE_EQ(results_host(ib), ubound - lbound);
    });
}

void batched_operator_2d()
{
    CoordBatch1 b1_min(0.0);
    CoordBatch1 b1_max(3.0);
    IVectBatch1 b1_ncells(4);
    CoordBatch2 b2_min(1.0);
    CoordBatch2 b2_max(2.0);
    IVectBatch2 b2_ncells(3);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IVectX x_ncells(16);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IVectY y_ncells(16);

    IDomainBatch1 gridb1 = ddc::init_discrete_space<IDimBatch1>(
            IDimBatch1::init<IDimBatch1>(b1_min, b1_max, b1_ncells));
    IDomainBatch2 gridb2 = ddc::init_discrete_space<IDimBatch2>(
            IDimBatch2::init<IDimBatch2>(b2_min, b2_max, b2_ncells));
    IDomainX gridx = ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_ncells));
    IDomainY gridy = ddc::init_discrete_space<IDimY>(IDimY::init<IDimY>(y_min, y_max, y_ncells));

    IDomainXY gridxy(gridx, gridy);

    DFieldXY quad_coeffs = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);

    Quadrature<IDomainXY, IDomainB1B2XY> quad_batched_operator(quad_coeffs.span_view());

    IDomainB1B2 gridb(gridb1, gridb2);

    DFieldB1B2 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            results.span_view(),
            KOKKOS_LAMBDA(IdxB1B2XY ibx) {
                double const x = ddc::coordinate(ddc::select<IDimX>(ibx));
                double const y = ddc::coordinate(ddc::select<IDimY>(ibx));
                double const b1 = ddc::coordinate(ddc::select<IDimBatch1>(ibx));
                double const b2 = ddc::coordinate(ddc::select<IDimBatch2>(ibx));
                return b1 * (x * y + 2 * x + b2 * y);
            });

    auto results_host = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), results.span_view());

    ddc::for_each(gridb, [&](IdxB1B2 ib) {
        double const b1 = ddc::coordinate(ddc::select<IDimBatch1>(ib));
        double const b2 = ddc::coordinate(ddc::select<IDimBatch2>(ib));
        double const exact
                = b1
                  * (y_max * y_max
                             * (2 * b2 * x_max - 2 * b2 * x_min + x_max * x_max - x_min * x_min)
                     + 4 * y_max * (x_max * x_max - x_min * x_min)
                     + y_min * y_min
                               * (-2 * b2 * x_max + 2 * b2 * x_min - x_max * x_max + x_min * x_min)
                     - 4 * y_min * (x_max * x_max - x_min * x_min))
                  / 4;
        EXPECT_NEAR(results_host(ib), exact, 1e-11);
    });
}

void batched_operator_1d_2d()
{
    CoordBatch1 b1_min(0.0);
    CoordBatch1 b1_max(3.0);
    IVectBatch1 b1_ncells(4);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IVectX x_ncells(16);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IVectY y_ncells(16);

    IDomainBatch1 gridb1 = ddc::init_discrete_space<IDimBatch1>(
            IDimBatch1::init<IDimBatch1>(b1_min, b1_max, b1_ncells));
    IDomainX gridx = ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_ncells));
    IDomainY gridy = ddc::init_discrete_space<IDimY>(IDimY::init<IDimY>(y_min, y_max, y_ncells));

    IDomainXY gridxy(gridx, gridy);

    DFieldXY quad_coeffs = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);

    Quadrature<IDomainXY, IDomainB1XY> quad_batched_operator(quad_coeffs.span_view());

    DFieldBatch1 results(gridb1);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            results.span_view(),
            KOKKOS_LAMBDA(IdxB1XY ibx) {
                double const x = ddc::coordinate(ddc::select<IDimX>(ibx));
                double const y = ddc::coordinate(ddc::select<IDimY>(ibx));
                double const b1 = ddc::coordinate(ddc::select<IDimBatch1>(ibx));
                return b1 * (x * y + 2 * x + 3 * y);
            });

    auto results_host = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), results.span_view());

    ddc::for_each(gridb1, [&](IdxBatch1 ib) {
        double const b1 = ddc::coordinate(ib);
        double const exact
                = b1
                  * (y_max * y_max * (x_max * x_max + 6 * x_max - x_min * x_min - 6 * x_min)
                     + 4 * y_max * (x_max * x_max - x_min * x_min)
                     + y_min * y_min * (-x_max * x_max - 6 * x_max + x_min * x_min + 6 * x_min)
                     - 4 * y_min * (x_max * x_max - x_min * x_min))
                  / 4;
        EXPECT_DOUBLE_EQ(results_host(ib), exact);
    });
}

void batched_operator_2d_1d()
{
    CoordBatch1 b1_min(0.0);
    CoordBatch1 b1_max(3.0);
    IVectBatch1 b1_ncells(4);
    CoordBatch2 b2_min(1.0);
    CoordBatch2 b2_max(2.0);
    IVectBatch2 b2_ncells(3);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IVectX x_ncells(16);

    IDomainBatch1 gridb1 = ddc::init_discrete_space<IDimBatch1>(
            IDimBatch1::init<IDimBatch1>(b1_min, b1_max, b1_ncells));
    IDomainBatch2 gridb2 = ddc::init_discrete_space<IDimBatch2>(
            IDimBatch2::init<IDimBatch2>(b2_min, b2_max, b2_ncells));
    IDomainX gridx = ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_ncells));

    DFieldX quad_coeffs = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx);

    Quadrature<IDomainX, IDomainB1B2X> quad_batched_operator(quad_coeffs.span_view());

    IDomainB1B2 gridb(gridb1, gridb2);

    DFieldB1B2 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            results.span_view(),
            KOKKOS_LAMBDA(IdxB1B2X ibx) {
                double const x = ddc::coordinate(ddc::select<IDimX>(ibx));
                double const b1 = ddc::coordinate(ddc::select<IDimBatch1>(ibx));
                double const b2 = ddc::coordinate(ddc::select<IDimBatch2>(ibx));
                return b1 * x + b2;
            });

    auto results_host = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), results.span_view());

    ddc::for_each(gridb, [&](IdxB1B2 ib) {
        double const b1 = ddc::coordinate(ddc::select<IDimBatch1>(ib));
        double const b2 = ddc::coordinate(ddc::select<IDimBatch2>(ib));
        double x = x_max;
        double const ubound = 0.5 * b1 * x * x + b2 * x;
        x = x_min;
        double const lbound = 0.5 * b1 * x * x + b2 * x;
        EXPECT_DOUBLE_EQ(results_host(ib), ubound - lbound);
    });
}

void batched_operator_2d_reordered()
{
    CoordBatch1 b1_min(0.0);
    CoordBatch1 b1_max(3.0);
    IVectBatch1 b1_ncells(4);
    CoordBatch2 b2_min(1.0);
    CoordBatch2 b2_max(2.0);
    IVectBatch2 b2_ncells(3);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IVectX x_ncells(16);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IVectY y_ncells(16);

    IDomainBatch1 gridb1 = ddc::init_discrete_space<IDimBatch1>(
            IDimBatch1::init<IDimBatch1>(b1_min, b1_max, b1_ncells));
    IDomainBatch2 gridb2 = ddc::init_discrete_space<IDimBatch2>(
            IDimBatch2::init<IDimBatch2>(b2_min, b2_max, b2_ncells));
    IDomainX gridx = ddc::init_discrete_space<IDimX>(IDimX::init<IDimX>(x_min, x_max, x_ncells));
    IDomainY gridy = ddc::init_discrete_space<IDimY>(IDimY::init<IDimY>(y_min, y_max, y_ncells));

    IDomainXY gridxy(gridx, gridy);

    DFieldXY quad_coeffs = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy);

    Quadrature<IDomainXY, IDomainXB1YB2> quad_batched_operator(quad_coeffs.span_view());

    IDomainB1B2 gridb(gridb1, gridb2);

    DFieldB1B2 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            results.span_view(),
            KOKKOS_LAMBDA(IdxXB1YB2 ibx) {
                double const x = ddc::coordinate(ddc::select<IDimX>(ibx));
                double const y = ddc::coordinate(ddc::select<IDimY>(ibx));
                double const b1 = ddc::coordinate(ddc::select<IDimBatch1>(ibx));
                double const b2 = ddc::coordinate(ddc::select<IDimBatch2>(ibx));
                return b1 * (x * y + 2 * x + b2 * y);
            });

    auto results_host = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), results.span_view());

    ddc::for_each(gridb, [&](IdxB1B2 ib) {
        double const b1 = ddc::coordinate(ddc::select<IDimBatch1>(ib));
        double const b2 = ddc::coordinate(ddc::select<IDimBatch2>(ib));
        double const exact
                = b1
                  * (y_max * y_max
                             * (2 * b2 * x_max - 2 * b2 * x_min + x_max * x_max - x_min * x_min)
                     + 4 * y_max * (x_max * x_max - x_min * x_min)
                     + y_min * y_min
                               * (-2 * b2 * x_max + 2 * b2 * x_min - x_max * x_max + x_min * x_min)
                     - 4 * y_min * (x_max * x_max - x_min * x_min))
                  / 4;
        EXPECT_NEAR(results_host(ib), exact, 1e-11);
    });
}

TEST(TrapezoidUniformNonPeriodicQuadrature, ExactForLinearBatch1D)
{
    batched_operator_1d();
}

TEST(TrapezoidUniformNonPeriodicQuadrature, ExactForLinearBatch1D2D)
{
    batched_operator_1d_2d();
}

TEST(TrapezoidUniformNonPeriodicQuadrature, ExactForLinearBatch2D1D)
{
    batched_operator_2d_1d();
}

TEST(TrapezoidUniformNonPeriodicQuadrature, ExactForLinearBatch2D)
{
    batched_operator_2d();
}

TEST(TrapezoidUniformNonPeriodicQuadrature, ExactForLinearBatch2DReordered)
{
    batched_operator_2d_reordered();
}
} // namespace
