// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

struct DDimX;
using DElemX = ddc::DiscreteElement<DDimX>;
using DVectX = ddc::DiscreteVector<DDimX>;
using DDomX = ddc::DiscreteDomain<DDimX>;

struct DDimY;
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

static DElemXY constexpr lbound_x_y {lbound_x, lbound_y};
static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);

TEST(ForEachSerialHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> view(storage.data(), dom);
    ddc::for_each(ddc::policies::serial_host, dom, [=](DElemX const ix) { view(ix) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachSerialHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> view(storage.data(), dom);
    ddc::for_each(ddc::policies::serial_host, dom, [=](DElemXY const ixy) { view(ixy) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachParallelHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> view(storage.data(), dom);
    ddc::for_each(ddc::policies::parallel_host, dom, [=](DElemX const ix) { view(ix) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachParallelHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> view(storage.data(), dom);
    ddc::for_each(ddc::policies::parallel_host, dom, [=](DElemXY const ixy) { view(ixy) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

static void TestForEachParallelDeviceOneDimension()
{
    DDomX const dom(lbound_x, nelems_x);
    ddc::Chunk<int, DDomX, ddc::DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ddc::ChunkSpan view(storage.span_view());
    ddc::for_each(
            ddc::policies::parallel_device,
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
    ddc::Chunk<int, DDomXY, ddc::DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ddc::ChunkSpan view(storage.span_view());
    ddc::for_each(
            ddc::policies::parallel_device,
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
