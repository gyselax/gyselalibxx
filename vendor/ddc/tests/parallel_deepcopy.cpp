// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <vector>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace DDC_HIP_5_7_ANONYMOUS_NAMESPACE_WORKAROUND(PARALLEL_DEEPCOPY_CPP)
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
    static DVectX constexpr nelems_x(2);

    static DElemY constexpr lbound_y(0);
    static DVectY constexpr nelems_y(2);

    static DElemXY constexpr lbound_x_y(lbound_x, lbound_y);
    static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);

} // namespace )

TEST(ParallelDeepcopy, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    ddc::Chunk chk(dom, ddc::HostAllocator<int>());
    chk(dom.front() + DVectXY(0, 0)) = 1;
    chk(dom.front() + DVectXY(1, 0)) = 2;
    chk(dom.front() + DVectXY(0, 1)) = 3;
    chk(dom.front() + DVectXY(1, 1)) = 4;
    ddc::Chunk chk_copy(dom, ddc::HostAllocator<int>());
    ddc::parallel_deepcopy(chk_copy, chk);
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(0, 0)), 1);
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(1, 0)), 2);
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(0, 1)), 3);
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(1, 1)), 4);
}

TEST(ParallelDeepcopy, TwoDimensionsWithExecutionSpace)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    ddc::Chunk chk(dom, ddc::HostAllocator<int>());
    chk(dom.front() + DVectXY(0, 0)) = 1;
    chk(dom.front() + DVectXY(1, 0)) = 2;
    chk(dom.front() + DVectXY(0, 1)) = 3;
    chk(dom.front() + DVectXY(1, 1)) = 4;
    ddc::Chunk chk_copy(dom, ddc::HostAllocator<int>());
    Kokkos::DefaultHostExecutionSpace const exec_space;
    ddc::parallel_deepcopy(exec_space, chk_copy, chk);
    exec_space.fence();
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(0, 0)), 1);
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(1, 0)), 2);
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(0, 1)), 3);
    EXPECT_EQ(chk_copy(dom.front() + DVectXY(1, 1)), 4);
}
