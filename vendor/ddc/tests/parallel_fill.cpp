// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <vector>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace DDC_HIP_5_7_ANONYMOUS_NAMESPACE_WORKAROUND(PARALLEL_FILL_CPP)
{
    struct DDimX
    {
    };
    using DElemX = ddc::DiscreteElement<DDimX>;
    using DVectX = ddc::DiscreteVector<DDimX>;
    using DDomX = ddc::DiscreteDomain<DDimX>;

    struct DDimY
    {
    };
    using DElemY = ddc::DiscreteElement<DDimY>;
    using DVectY = ddc::DiscreteVector<DDimY>;
    using DDomY = ddc::DiscreteDomain<DDimY>;

    using DElemXY = ddc::DiscreteElement<DDimX, DDimY>;
    using DVectXY = ddc::DiscreteVector<DDimX, DDimY>;
    using DDomXY = ddc::DiscreteDomain<DDimX, DDimY>;

    static DElemX constexpr lbound_x(0);
    static DVectX constexpr nelems_x(10);

    static DElemY constexpr lbound_y(0);
    static DVectY constexpr nelems_y(12);

    static DElemXY constexpr lbound_x_y(lbound_x, lbound_y);
    static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);

} // namespace )

TEST(ParallelFill, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> const view(storage.data(), dom);
    ddc::parallel_fill(view, 1);
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ParallelFill, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> const view(storage.data(), dom);
    ddc::parallel_fill(view, 1);
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ParallelFill, OneDimensionWithExecutionSpace)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> const view(storage.data(), dom);
    Kokkos::DefaultHostExecutionSpace const exec_space;
    ddc::parallel_fill(exec_space, view, 1);
    exec_space.fence();
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ParallelFill, TwoDimensionsWithExecutionSpace)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> const view(storage.data(), dom);
    Kokkos::DefaultHostExecutionSpace const exec_space;
    ddc::parallel_fill(exec_space, view, 1);
    exec_space.fence();
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}
