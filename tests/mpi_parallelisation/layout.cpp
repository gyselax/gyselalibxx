// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "mpilayout.hpp"

namespace {

struct X
{
};
struct Y
{
};

struct GridX : UniformGridBase<X>
{
};
struct GridY : UniformGridBase<Y>
{
};

using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;

using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;
using IdxStepXY = IdxStep<GridX, GridY>;

using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeYX = IdxRange<GridY, GridX>;

using CoordX = Coord<X>;
using CoordY = Coord<Y>;

using XYDistribLayout = MPILayout<IdxRangeXY, GridX, GridY>;
using YDistribLayout = MPILayout<IdxRangeXY, GridY>;

TEST(Layout, MinimalDomainDistribution)
{
    IdxStepX const x_size(8);
    IdxStepY const y_size(16);

    IdxXY const idx_range_start(0, 0);
    IdxStepXY const idx_range_size(x_size, y_size);
    IdxRangeXY const global_idx_range(idx_range_start, idx_range_size);

    int n_procs = 4;
    int expected_local_x_extent = x_size.value() / n_procs;

    IdxX const idx_range_x_start(idx_range_start);

    XYDistribLayout layout;
    for (int i(0); i < n_procs; ++i) {
        IdxRangeXY const local_idx_range
                = layout.distribute_idx_range(global_idx_range, n_procs, i);
        EXPECT_EQ(local_idx_range.extent<GridX>().value(), expected_local_x_extent);
        EXPECT_EQ(
                ddc::select<GridX>(local_idx_range.front()),
                idx_range_x_start + i * expected_local_x_extent);
        EXPECT_EQ(local_idx_range.extent<GridY>().value(), y_size.value());
    }
}

TEST(Layout, SpreadDomainDistribution)
{
    IdxStepX const x_size(10);
    IdxStepY const y_size(15);

    IdxXY const idx_range_start(0, 0);
    IdxStepXY const idx_range_size(x_size, y_size);
    IdxRangeXY const global_idx_range(idx_range_start, idx_range_size);

    int const n_procs = 6;
    int const expected_procs_x = 2;
    int const expected_procs_y = 3;
    int const expected_local_x_extent = x_size / expected_procs_x;
    int const expected_local_y_extent = y_size / expected_procs_y;

    IdxX const idx_range_x_start(idx_range_start);
    IdxY const idx_range_y_start(idx_range_start);

    XYDistribLayout layout;
    for (int i(0); i < n_procs; ++i) {
        IdxRangeXY const local_idx_range
                = layout.distribute_idx_range(global_idx_range, n_procs, i);
        EXPECT_EQ(local_idx_range.extent<GridX>().value(), expected_local_x_extent);
        EXPECT_EQ(local_idx_range.extent<GridY>().value(), expected_local_y_extent);
        int const x_rank = i / expected_procs_y;
        int const y_rank = i % expected_procs_y;
        EXPECT_EQ(
                ddc::select<GridX>(local_idx_range.front()),
                idx_range_x_start + x_rank * expected_local_x_extent);
        EXPECT_EQ(
                ddc::select<GridY>(local_idx_range.front()),
                idx_range_y_start + y_rank * expected_local_y_extent);
    }
}

TEST(Layout, DomainSelection)
{
    IdxStepX const x_size(10);
    IdxStepY const y_size(15);

    IdxXY const idx_range_start(0, 0);
    IdxStepXY const idx_range_size(x_size, y_size);
    IdxRangeXY const global_idx_range(idx_range_start, idx_range_size);

    int n_procs = 5;
    int expected_local_y_extent = 3;

    IdxY const idx_range_y_start(idx_range_start);

    YDistribLayout layout;
    for (int i(0); i < n_procs; ++i) {
        IdxRangeXY const local_idx_range
                = layout.distribute_idx_range(global_idx_range, n_procs, i);
        EXPECT_EQ(local_idx_range.extent<GridX>(), global_idx_range.extent<GridX>());
        EXPECT_EQ(
                ddc::select<GridY>(local_idx_range.front()),
                idx_range_y_start + i * expected_local_y_extent);
        EXPECT_EQ(local_idx_range.extent<GridY>().value(), expected_local_y_extent);
    }
}
} // namespace
