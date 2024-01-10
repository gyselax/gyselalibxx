// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

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

TEST(TransformReduceSerialHost, ZeroDimension)
{
    DDom0D const dom;
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDom0D> chunk(storage.data(), dom);
    int count = 0;
    ddc::for_each(dom, [&](DElem0D const i) { chunk(i) = count++; });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::serial_host,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    [&](DElem0D const i) { return chunk(i); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceSerialHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> chunk(storage.data(), dom);
    int count = 0;
    ddc::for_each(dom, [&](DElemX const ix) { chunk(ix) = count++; });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::serial_host,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    [&](DElemX const ix) { return chunk(ix); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceSerialHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> chunk(storage.data(), dom);
    int count = 0;
    ddc::for_each(dom, [&](DElemXY const ixy) { chunk(ixy) = count++; });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::serial_host,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    [&](DElemXY const ixy) { return chunk(ixy); }),
            dom.size() * (dom.size() - 1) / 2);
}
TEST(TransformReduceParallelHost, ZeroDimension)
{
    DDom0D const dom;
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDom0D> chunk(storage.data(), dom);
    int count = 0;
    ddc::for_each(dom, [&](DElem0D const i) { chunk(i) = count++; });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::parallel_host,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    [&](DElem0D const i) { return chunk(i); }),
            dom.size() * (dom.size() - 1) / 2);
}
TEST(TransformReduceParallelHost, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomX> chunk(storage.data(), dom);
    int count = 0;
    ddc::for_each(dom, [&](DElemX const ix) { chunk(ix) = count++; });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::parallel_host,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    [&](DElemX const ix) { return chunk(ix); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceParallelHost, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ddc::ChunkSpan<int, DDomXY> chunk(storage.data(), dom);
    int count = 0;
    ddc::for_each(dom, [&](DElemXY const ixy) { chunk(ixy) = count++; });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::parallel_host,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    [&](DElemXY const ixy) { return chunk(ixy); }),
            dom.size() * (dom.size() - 1) / 2);
}

static void TestTransformReduceParallelDeviceZeroDimension()
{
    DDom0D const dom;
    ddc::Chunk<int, DDom0D, ddc::DeviceAllocator<int>> storage(dom);
    ddc::ChunkSpan const chunk(storage.span_view());
    Kokkos::View<int> count("count");
    Kokkos::deep_copy(count, 0);
    ddc::for_each(
            ddc::policies::parallel_device,
            dom,
            KOKKOS_LAMBDA(DElem0D const i) { chunk(i) = Kokkos::atomic_fetch_add(&count(), 1); });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::parallel_device,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    KOKKOS_LAMBDA(DElem0D const i) { return chunk(i); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceParallelDevice, ZeroDimension)
{
    TestTransformReduceParallelDeviceZeroDimension();
}

static void TestTransformReduceParallelDeviceOneDimension()
{
    DDomX const dom(lbound_x, nelems_x);
    ddc::Chunk<int, DDomX, ddc::DeviceAllocator<int>> storage(dom);
    ddc::ChunkSpan const chunk(storage.span_view());
    Kokkos::View<int> count("count");
    Kokkos::deep_copy(count, 0);
    ddc::for_each(
            ddc::policies::parallel_device,
            dom,
            KOKKOS_LAMBDA(DElemX const ix) { chunk(ix) = Kokkos::atomic_fetch_add(&count(), 1); });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::parallel_device,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    KOKKOS_LAMBDA(DElemX const ix) { return chunk(ix); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceParallelDevice, OneDimension)
{
    TestTransformReduceParallelDeviceOneDimension();
}

static void TestTransformReduceParallelDeviceTwoDimensions()
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    ddc::Chunk<int, DDomXY, ddc::DeviceAllocator<int>> storage(dom);
    ddc::ChunkSpan const chunk(storage.span_view());
    Kokkos::View<int> count("count");
    Kokkos::deep_copy(count, 0);
    ddc::for_each(
            ddc::policies::parallel_device,
            dom,
            KOKKOS_LAMBDA(DElemXY const ixy) {
                chunk(ixy) = Kokkos::atomic_fetch_add(&count(), 1);
            });
    EXPECT_EQ(
            ddc::transform_reduce(
                    ddc::policies::parallel_device,
                    dom,
                    0,
                    ddc::reducer::sum<int>(),
                    KOKKOS_LAMBDA(DElemXY const ixy) { return chunk(ixy); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceParallelDevice, TwoDimensions)
{
    TestTransformReduceParallelDeviceTwoDimensions();
}
