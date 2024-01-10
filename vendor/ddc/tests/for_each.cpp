// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

using DElem0D = ddc::DiscreteElement<>;
using DVect0D = ddc::DiscreteVector<>;
using DDom0D = ddc::DiscreteDomain<>;

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

static DElemXY constexpr lbound_x_y(lbound_x, lbound_y);
static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);

} // namespace

TEST(ForEachSerialHost, Empty)
{
    DDomX const dom(lbound_x, DVectX(0));
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> view(storage.data(), dom);
    ddc::for_each(ddc::policies::serial_host, dom, [=](DElemX const ix) { view(ix) += 1; });
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size())
            << std::count(storage.begin(), storage.end(), 1) << std::endl;
}

TEST(ForEachSerialHost, ZeroDimension)
{
    DDom0D const dom;
    int storage = 0;
    ddc::ChunkSpan<int, DDom0D> view(&storage, dom);
    ddc::for_each(ddc::policies::serial_host, dom, [=](DElem0D const ii) { view(ii) += 1; });
    EXPECT_EQ(storage, 1) << storage << std::endl;
}

TEST(ForEachSerialHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> view(storage.data(), dom);
    ddc::for_each(ddc::policies::serial_host, dom, [=](DElemX const ix) { view(ix) += 1; });
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachSerialHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> view(storage.data(), dom);
    ddc::for_each(ddc::policies::serial_host, dom, [=](DElemXY const ixy) { view(ixy) += 1; });
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachParallelHost, ZeroDimension)
{
    DDom0D const dom;
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDom0D> view(storage.data(), dom);
    ddc::for_each(ddc::policies::parallel_host, dom, [=](DElem0D const i) { view(i) += 1; });
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachParallelHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> view(storage.data(), dom);
    ddc::for_each(ddc::policies::parallel_host, dom, [=](DElemX const ix) { view(ix) += 1; });
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachParallelHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> view(storage.data(), dom);
    ddc::for_each(ddc::policies::parallel_host, dom, [=](DElemXY const ixy) { view(ixy) += 1; });
    EXPECT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

namespace {

void TestForEachParallelDeviceZeroDimension()
{
    DDom0D const dom;
    ddc::Chunk<int, DDom0D, ddc::DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ddc::ChunkSpan view(storage.span_view());
    ddc::for_each(
            ddc::policies::parallel_device,
            dom,
            KOKKOS_LAMBDA(DElem0D const i) { view(i) += 1; });
    int const* const ptr = storage.data_handle();
    int sum;
    Kokkos::parallel_reduce(
            dom.size(),
            KOKKOS_LAMBDA(std::size_t i, int& local_sum) { local_sum += ptr[i]; },
            Kokkos::Sum<int>(sum));
    EXPECT_EQ(sum, dom.size());
}

} // namespace

TEST(ForEachParallelDevice, ZeroDimension)
{
    TestForEachParallelDeviceZeroDimension();
}

namespace {

void TestForEachParallelDeviceOneDimension()
{
    DDomX const dom(lbound_x, nelems_x);
    ddc::Chunk<int, DDomX, ddc::DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ddc::ChunkSpan view(storage.span_view());
    ddc::for_each(
            ddc::policies::parallel_device,
            dom,
            KOKKOS_LAMBDA(DElemX const ix) { view(ix) += 1; });
    int const* const ptr = storage.data_handle();
    int sum;
    Kokkos::parallel_reduce(
            dom.size(),
            KOKKOS_LAMBDA(std::size_t i, int& local_sum) { local_sum += ptr[i]; },
            Kokkos::Sum<int>(sum));
    EXPECT_EQ(sum, dom.size());
}

} // namespace

TEST(ForEachParallelDevice, OneDimension)
{
    TestForEachParallelDeviceOneDimension();
}

namespace {

void TestForEachParallelDeviceTwoDimensions()
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    ddc::Chunk<int, DDomXY, ddc::DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ddc::ChunkSpan view(storage.span_view());
    ddc::for_each(
            ddc::policies::parallel_device,
            dom,
            KOKKOS_LAMBDA(DElemXY const ixy) { view(ixy) += 1; });
    int const* const ptr = storage.data_handle();
    int sum;
    Kokkos::parallel_reduce(
            dom.size(),
            KOKKOS_LAMBDA(std::size_t i, int& local_sum) { local_sum += ptr[i]; },
            Kokkos::Sum<int>(sum));
    EXPECT_EQ(sum, dom.size());
}

} // namespace

TEST(ForEachParallelDevice, TwoDimensions)
{
    TestForEachParallelDeviceTwoDimensions();
}
