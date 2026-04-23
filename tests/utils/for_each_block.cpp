#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <ddc_for_each_block.hpp>

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

DElemX constexpr lbound_x = ddc::init_trivial_half_bounded_space<DDimX>();
DVectX constexpr nelems_x(10);

DElemY constexpr lbound_y = ddc::init_trivial_half_bounded_space<DDimY>();
DVectY constexpr nelems_y(12);

DElemXY constexpr lbound_x_y(lbound_x, lbound_y);
DVectXY constexpr nelems_x_y(nelems_x, nelems_y);

TEST(ForEachSerialBlockHost, IsPowerOf2)
{
    EXPECT_FALSE(ddcHelper::detail::is_power_of_2(0));
    EXPECT_TRUE(ddcHelper::detail::is_power_of_2(1));
    EXPECT_TRUE(ddcHelper::detail::is_power_of_2(2));
    EXPECT_FALSE(ddcHelper::detail::is_power_of_2(3));
    EXPECT_TRUE(ddcHelper::detail::is_power_of_2(4));
}

TEST(DistributeBlocks, D)
{
    {
        std::array<ddc::DiscreteVectorElement, 1> size {10};
        std::array<ddc::DiscreteVectorElement, 1> nb_blocks_per_dim;
        ddcHelper::detail::distribute_blocks(2, size, nb_blocks_per_dim);
        EXPECT_EQ(nb_blocks_per_dim, (std::array<ddc::DiscreteVectorElement, 1> {2}));
    }
    {
        std::array<ddc::DiscreteVectorElement, 1> size {10};
        std::array<ddc::DiscreteVectorElement, 1> nb_blocks_per_dim;
        ddcHelper::detail::distribute_blocks(8, size, nb_blocks_per_dim);
        EXPECT_EQ(nb_blocks_per_dim, (std::array<ddc::DiscreteVectorElement, 1> {8}));
    }
    {
        std::array<ddc::DiscreteVectorElement, 3> size {3, 4, 5};
        std::array<ddc::DiscreteVectorElement, 3> nb_blocks_per_dim;
        ddcHelper::detail::distribute_blocks(32, size, nb_blocks_per_dim);
        EXPECT_EQ(nb_blocks_per_dim, (std::array<ddc::DiscreteVectorElement, 3> {2, 4, 4}));
    }
}

TEST(ForEachSerialBlockHost, OneDimension)
{
    for (ddc::DiscreteVectorElement const nb_blocks : {1, 2, 4, 8}) {
        DDomX const dom(lbound_x, nelems_x);
        ddc::Chunk elems_count("count", dom, ddc::KokkosAllocator<int, Kokkos::HostSpace>());
        int measured_nb_blocks = 0;
        ddc::parallel_fill(elems_count, 0);
        ddcHelper::host_for_each_block(dom, nb_blocks, [&](DDomX const domx) {
            ddc::host_for_each(domx, [&](DElemX ix) { elems_count(ix)++; });
            ++measured_nb_blocks;
        });
        ddc::host_for_each(dom, [&](DElemX ix) { EXPECT_EQ(elems_count(ix), 1); });
        EXPECT_EQ(measured_nb_blocks, nb_blocks);
    }
}

TEST(ForEachSerialBlockHost, TwoDimensions)
{
    for (ddc::DiscreteVectorElement const nb_blocks : {1, 2, 4, 8}) {
        DDomXY const dom(lbound_x_y, nelems_x_y);
        ddc::Chunk elems_count("count", dom, ddc::KokkosAllocator<int, Kokkos::HostSpace>());
        int measured_nb_blocks = 0;
        ddc::parallel_fill(elems_count, 0);
        ddcHelper::host_for_each_block(dom, nb_blocks, [&](DDomXY const domxy) {
            ddc::host_for_each(domxy, [&](DElemXY ixy) { elems_count(ixy)++; });
            ++measured_nb_blocks;
        });
        ddc::host_for_each(dom, [&](DElemXY ixy) { EXPECT_EQ(elems_count(ixy), 1); });
        EXPECT_EQ(measured_nb_blocks, nb_blocks);
    }
}
