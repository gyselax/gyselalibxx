// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

struct DDimX;
using ElemX = DiscreteCoordinate<DDimX>;
using DVectX = DiscreteVector<DDimX>;
using DDomX = DiscreteDomain<DDimX>;

struct DDimY;
using ElemY = DiscreteCoordinate<DDimY>;
using DVectY = DiscreteVector<DDimY>;
using DDomY = DiscreteDomain<DDimY>;

using ElemXY = DiscreteCoordinate<DDimX, DDimY>;
using DVectXY = DiscreteVector<DDimX, DDimY>;
using DDomXY = DiscreteDomain<DDimX, DDimY>;

static ElemX constexpr lbound_x(0);
static DVectX constexpr nelems_x(10);

static ElemY constexpr lbound_y(0);
static DVectY constexpr nelems_y(12);

static ElemXY constexpr lbound_x_y {lbound_x, lbound_y};
static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);

TEST(TransformReduceSerial, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomX> chunk(storage.data(), dom);
    int count = 0;
    for_each(dom, [&](ElemX const ix) { chunk(ix) = count++; });
    ASSERT_EQ(
            transform_reduce(
                    policies::serial,
                    dom,
                    0,
                    reducer::sum<int>(),
                    [&](ElemX const ix) { return chunk(ix); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceSerial, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomXY> chunk(storage.data(), dom);
    int count = 0;
    for_each(dom, [&](ElemXY const ixy) { chunk(ixy) = count++; });
    ASSERT_EQ(
            transform_reduce(
                    policies::serial,
                    dom,
                    0,
                    reducer::sum<int>(),
                    [&](ElemXY const ixy) { return chunk(ixy); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceOmp, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomX> chunk(storage.data(), dom);
    int count = 0;
    for_each(dom, [&](ElemX const ix) { chunk(ix) = count++; });
    ASSERT_EQ(
            transform_reduce(
                    policies::omp,
                    dom,
                    0,
                    reducer::sum<int>(),
                    [&](ElemX const ix) { return chunk(ix); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceOmp, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomXY> chunk(storage.data(), dom);
    int count = 0;
    for_each(dom, [&](ElemXY const ixy) { chunk(ixy) = count++; });
    ASSERT_EQ(
            transform_reduce(
                    policies::omp,
                    dom,
                    0,
                    reducer::sum<int>(),
                    [&](ElemXY const ixy) { return chunk(ixy); }),
            dom.size() * (dom.size() - 1) / 2);
}

static void TestTransformReduceKokkosOneDimension()
{
    DDomX const dom(lbound_x, nelems_x);
    Chunk<int, DDomX, DeviceAllocator<int>> storage(dom);
    ChunkSpan const chunk(storage.span_view());
    Kokkos::View<int> count("count");
    Kokkos::deep_copy(count, 0);
    for_each(
            policies::parallel_device,
            dom,
            DDC_LAMBDA(ElemX const ix) { chunk(ix) = Kokkos::atomic_fetch_add(&count(), 1); });
    ASSERT_EQ(
            transform_reduce(
                    policies::parallel_device,
                    dom,
                    0,
                    reducer::sum<int>(),
                    DDC_LAMBDA(ElemX const ix) { return chunk(ix); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceKokkos, OneDimension)
{
    TestTransformReduceKokkosOneDimension();
}

static void TestTransformReduceKokkosTwoDimensions()
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    Chunk<int, DDomXY, DeviceAllocator<int>> storage(dom);
    ChunkSpan const chunk(storage.span_view());
    Kokkos::View<int> count("count");
    Kokkos::deep_copy(count, 0);
    for_each(
            policies::parallel_device,
            dom,
            DDC_LAMBDA(ElemXY const ixy) { chunk(ixy) = Kokkos::atomic_fetch_add(&count(), 1); });
    ASSERT_EQ(
            transform_reduce(
                    policies::parallel_device,
                    dom,
                    0,
                    reducer::sum<int>(),
                    DDC_LAMBDA(ElemXY const ixy) { return chunk(ixy); }),
            dom.size() * (dom.size() - 1) / 2);
}

TEST(TransformReduceKokkos, TwoDimensions)
{
    TestTransformReduceKokkosTwoDimensions();
}
