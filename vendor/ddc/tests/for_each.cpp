// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

struct DDimX;
using DElemX = DiscreteElement<DDimX>;
using DVectX = DiscreteVector<DDimX>;
using DDomX = DiscreteDomain<DDimX>;

struct DDimY;
using DElemY = DiscreteElement<DDimY>;
using DVectY = DiscreteVector<DDimY>;
using DDomY = DiscreteDomain<DDimY>;

using DElemXY = DiscreteElement<DDimX, DDimY>;
using DVectXY = DiscreteVector<DDimX, DDimY>;
using DDomXY = DiscreteDomain<DDimX, DDimY>;

static DElemX constexpr lbound_x(0);
static DVectX constexpr nelems_x(10);

static DElemY constexpr lbound_y(0);
static DVectY constexpr nelems_y(12);

static DElemXY constexpr lbound_x_y {lbound_x, lbound_y};
static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);

TEST(ForEachSerialHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomX> view(storage.data(), dom);
    for_each(policies::serial_host, dom, [=](DElemX const ix) { view(ix) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachSerialHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomXY> view(storage.data(), dom);
    for_each(policies::serial_host, dom, [=](DElemXY const ixy) { view(ixy) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachParallelHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomX> view(storage.data(), dom);
    for_each(policies::parallel_host, dom, [=](DElemX const ix) { view(ix) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachParallelHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomXY> view(storage.data(), dom);
    for_each(policies::parallel_host, dom, [=](DElemXY const ixy) { view(ixy) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

static void TestForEachParallelDeviceOneDimension()
{
    DDomX const dom(lbound_x, nelems_x);
    Chunk<int, DDomX, DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ChunkSpan view(storage.span_view());
    for_each(
            policies::parallel_device,
            dom,
            DDC_LAMBDA(DElemX const ix) { view(ix) += 1; });
    int const* const ptr = storage.data();
    int sum;
    Kokkos::parallel_reduce(
            dom.size(),
            KOKKOS_LAMBDA(std::size_t i, int& local_sum) { local_sum += ptr[i]; },
            Kokkos::Sum<int>(sum));
    ASSERT_EQ(sum, dom.size());
}

TEST(ForEachParallelDevice, OneDimension)
{
    TestForEachParallelDeviceOneDimension();
}

static void TestForEachParallelDeviceTwoDimensions()
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    Chunk<int, DDomXY, DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ChunkSpan view(storage.span_view());
    for_each(
            policies::parallel_device,
            dom,
            DDC_LAMBDA(DElemXY const ixy) { view(ixy) += 1; });
    int const* const ptr = storage.data();
    int sum;
    Kokkos::parallel_reduce(
            dom.size(),
            KOKKOS_LAMBDA(std::size_t i, int& local_sum) { local_sum += ptr[i]; },
            Kokkos::Sum<int>(sum));
    ASSERT_EQ(sum, dom.size());
}

TEST(ForEachParallelDevice, TwoDimensions)
{
    TestForEachParallelDeviceTwoDimensions();
}
