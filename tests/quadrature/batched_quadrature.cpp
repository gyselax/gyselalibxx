// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
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

using CoordBatch1 = const Coord<Batch1>;
using CoordBatch2 = const Coord<Batch2>;
using CoordX = const Coord<X>;
using CoordY = const Coord<Y>;

struct GridBatch1 : UniformGridBase<Batch1>
{
};
struct GridBatch2 : UniformGridBase<Batch2>
{
};

struct GridX : UniformGridBase<X>
{
};
struct GridY : UniformGridBase<Y>
{
};

using IdxBatch1 = Idx<GridBatch1>;
using IdxBatch2 = Idx<GridBatch2>;
using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;
using IdxB1X = Idx<GridBatch1, GridX>;
using IdxB1B2 = Idx<GridBatch1, GridBatch2>;
using IdxB1XY = Idx<GridBatch1, GridX, GridY>;
using IdxB1B2X = Idx<GridBatch1, GridBatch2, GridX>;
using IdxB1B2XY = Idx<GridBatch1, GridBatch2, GridX, GridY>;
using IdxXB1YB2 = Idx<GridX, GridBatch1, GridY, GridBatch2>;

using IdxStepBatch1 = const IdxStep<GridBatch1>;
using IdxStepBatch2 = const IdxStep<GridBatch2>;
using IdxStepX = const IdxStep<GridX>;
using IdxStepY = const IdxStep<GridY>;

using IdxRangeBatch1 = IdxRange<GridBatch1>;
using IdxRangeBatch2 = IdxRange<GridBatch2>;
using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeB1X = IdxRange<GridBatch1, GridX>;
using IdxRangeB1XY = IdxRange<GridBatch1, GridX, GridY>;
using IdxRangeB1B2X = IdxRange<GridBatch1, GridBatch2, GridX>;
using IdxRangeB1B2XY = IdxRange<GridBatch1, GridBatch2, GridX, GridY>;
using IdxRangeXB1YB2 = IdxRange<GridX, GridBatch1, GridY, GridBatch2>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeB1B2 = IdxRange<GridBatch1, GridBatch2>;

using DFieldMemBatch1 = DFieldMem<IdxRangeBatch1>;
using DFieldMemX = DFieldMem<IdxRangeX>;
using DFieldMemY = DFieldMem<IdxRangeY>;
using DFieldMemXY = DFieldMem<IdxRangeXY>;
using DFieldMemB1B2 = DFieldMem<IdxRangeB1B2>;

void batched_operator_1d()
{
    CoordBatch1 b_min(0.0);
    CoordBatch1 b_max(3.0);
    IdxStepBatch1 b_ncells(4);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IdxStepX x_ncells(16);

    IdxRangeBatch1 gridb = ddc::init_discrete_space<GridBatch1>(
            GridBatch1::init<GridBatch1>(b_min, b_max, b_ncells));
    IdxRangeX gridx = ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_ncells));

    DFieldMemX quad_coeffs(trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx));

    Quadrature<IdxRangeX, IdxRangeB1X> quad_batched_operator(get_field(quad_coeffs));

    DFieldMemBatch1 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            get_field(results),
            KOKKOS_LAMBDA(IdxB1X ibx) {
                double b = ddc::coordinate(ddc::select<GridBatch1>(ibx));
                double x = ddc::coordinate(ddc::select<GridX>(ibx));
                return b * x + 2;
            });

    auto results_host = ddc::create_mirror_view_and_copy(get_field(results));

    ddc::for_each(gridb, [&](IdxBatch1 ib) {
        double b = ddc::coordinate(ddc::select<GridBatch1>(ib));
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
    IdxStepBatch1 b1_ncells(4);
    CoordBatch2 b2_min(1.0);
    CoordBatch2 b2_max(2.0);
    IdxStepBatch2 b2_ncells(3);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IdxStepX x_ncells(16);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IdxStepY y_ncells(16);

    IdxRangeBatch1 gridb1 = ddc::init_discrete_space<GridBatch1>(
            GridBatch1::init<GridBatch1>(b1_min, b1_max, b1_ncells));
    IdxRangeBatch2 gridb2 = ddc::init_discrete_space<GridBatch2>(
            GridBatch2::init<GridBatch2>(b2_min, b2_max, b2_ncells));
    IdxRangeX gridx = ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_ncells));
    IdxRangeY gridy = ddc::init_discrete_space<GridY>(GridY::init<GridY>(y_min, y_max, y_ncells));

    IdxRangeXY gridxy(gridx, gridy);

    DFieldMemXY quad_coeffs(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy));
    Quadrature<IdxRangeXY, IdxRangeB1B2XY> quad_batched_operator(get_field(quad_coeffs));

    IdxRangeB1B2 gridb(gridb1, gridb2);

    DFieldMemB1B2 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            get_field(results),
            KOKKOS_LAMBDA(IdxB1B2XY ibx) {
                double const x = ddc::coordinate(ddc::select<GridX>(ibx));
                double const y = ddc::coordinate(ddc::select<GridY>(ibx));
                double const b1 = ddc::coordinate(ddc::select<GridBatch1>(ibx));
                double const b2 = ddc::coordinate(ddc::select<GridBatch2>(ibx));
                return b1 * (x * y + 2 * x + b2 * y);
            });

    auto results_host = ddc::create_mirror_view_and_copy(get_field(results));

    ddc::for_each(gridb, [&](IdxB1B2 ib) {
        double const b1 = ddc::coordinate(ddc::select<GridBatch1>(ib));
        double const b2 = ddc::coordinate(ddc::select<GridBatch2>(ib));
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
    IdxStepBatch1 b1_ncells(4);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IdxStepX x_ncells(16);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IdxStepY y_ncells(16);

    IdxRangeBatch1 gridb1 = ddc::init_discrete_space<GridBatch1>(
            GridBatch1::init<GridBatch1>(b1_min, b1_max, b1_ncells));
    IdxRangeX gridx = ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_ncells));
    IdxRangeY gridy = ddc::init_discrete_space<GridY>(GridY::init<GridY>(y_min, y_max, y_ncells));

    IdxRangeXY gridxy(gridx, gridy);

    DFieldMemXY quad_coeffs(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy));

    Quadrature<IdxRangeXY, IdxRangeB1XY> quad_batched_operator(get_field(quad_coeffs));

    DFieldMemBatch1 results(gridb1);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            get_field(results),
            KOKKOS_LAMBDA(IdxB1XY ibx) {
                double const x = ddc::coordinate(ddc::select<GridX>(ibx));
                double const y = ddc::coordinate(ddc::select<GridY>(ibx));
                double const b1 = ddc::coordinate(ddc::select<GridBatch1>(ibx));
                return b1 * (x * y + 2 * x + 3 * y);
            });

    auto results_host = ddc::create_mirror_view_and_copy(get_field(results));

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
    IdxStepBatch1 b1_ncells(4);
    CoordBatch2 b2_min(1.0);
    CoordBatch2 b2_max(2.0);
    IdxStepBatch2 b2_ncells(3);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IdxStepX x_ncells(16);

    IdxRangeBatch1 gridb1 = ddc::init_discrete_space<GridBatch1>(
            GridBatch1::init<GridBatch1>(b1_min, b1_max, b1_ncells));
    IdxRangeBatch2 gridb2 = ddc::init_discrete_space<GridBatch2>(
            GridBatch2::init<GridBatch2>(b2_min, b2_max, b2_ncells));
    IdxRangeX gridx = ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_ncells));

    DFieldMemX quad_coeffs(trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx));
    Quadrature<IdxRangeX, IdxRangeB1B2X> quad_batched_operator(get_field(quad_coeffs));

    IdxRangeB1B2 gridb(gridb1, gridb2);

    DFieldMemB1B2 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            get_field(results),
            KOKKOS_LAMBDA(IdxB1B2X ibx) {
                double const x = ddc::coordinate(ddc::select<GridX>(ibx));
                double const b1 = ddc::coordinate(ddc::select<GridBatch1>(ibx));
                double const b2 = ddc::coordinate(ddc::select<GridBatch2>(ibx));
                return b1 * x + b2;
            });

    auto results_host = ddc::create_mirror_view_and_copy(get_field(results));

    ddc::for_each(gridb, [&](IdxB1B2 ib) {
        double const b1 = ddc::coordinate(ddc::select<GridBatch1>(ib));
        double const b2 = ddc::coordinate(ddc::select<GridBatch2>(ib));
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
    IdxStepBatch1 b1_ncells(4);
    CoordBatch2 b2_min(1.0);
    CoordBatch2 b2_max(2.0);
    IdxStepBatch2 b2_ncells(3);
    CoordX x_min(4.0);
    CoordX x_max(8.0);
    IdxStepX x_ncells(16);
    CoordY y_min(4.0);
    CoordY y_max(8.0);
    IdxStepY y_ncells(16);

    IdxRangeBatch1 gridb1 = ddc::init_discrete_space<GridBatch1>(
            GridBatch1::init<GridBatch1>(b1_min, b1_max, b1_ncells));
    IdxRangeBatch2 gridb2 = ddc::init_discrete_space<GridBatch2>(
            GridBatch2::init<GridBatch2>(b2_min, b2_max, b2_ncells));
    IdxRangeX gridx = ddc::init_discrete_space<GridX>(GridX::init<GridX>(x_min, x_max, x_ncells));
    IdxRangeY gridy = ddc::init_discrete_space<GridY>(GridY::init<GridY>(y_min, y_max, y_ncells));

    IdxRangeXY gridxy(gridx, gridy);

    DFieldMemXY quad_coeffs(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridxy));
    Quadrature<IdxRangeXY, IdxRangeXB1YB2> quad_batched_operator(get_field(quad_coeffs));

    IdxRangeB1B2 gridb(gridb1, gridb2);

    DFieldMemB1B2 results(gridb);
    quad_batched_operator(
            Kokkos::DefaultExecutionSpace(),
            get_field(results),
            KOKKOS_LAMBDA(IdxXB1YB2 ibx) {
                double const x = ddc::coordinate(ddc::select<GridX>(ibx));
                double const y = ddc::coordinate(ddc::select<GridY>(ibx));
                double const b1 = ddc::coordinate(ddc::select<GridBatch1>(ibx));
                double const b2 = ddc::coordinate(ddc::select<GridBatch2>(ibx));
                return b1 * (x * y + 2 * x + b2 * y);
            });

    auto results_host = ddc::create_mirror_view_and_copy(get_field(results));

    ddc::for_each(gridb, [&](IdxB1B2 ib) {
        double const b1 = ddc::coordinate(ddc::select<GridBatch1>(ib));
        double const b2 = ddc::coordinate(ddc::select<GridBatch2>(ib));
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
