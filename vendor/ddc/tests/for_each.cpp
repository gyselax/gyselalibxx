// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

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

TEST(ForEachSerial, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomX> view(storage.data(), dom);
    for_each(policies::serial, dom, [=](ElemX const ix) { view(ix) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachSerial, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomXY> view(storage.data(), dom);
    for_each(policies::serial, dom, [=](ElemXY const ixy) { view(ixy) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachOmp, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomX> view(storage.data(), dom);
    for_each(policies::omp, dom, [=](ElemX const ix) { view(ix) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

TEST(ForEachOmp, TwoDimensions)
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    std::vector<int> storage(dom.size(), 0);
    ChunkSpan<int, DDomXY> view(storage.data(), dom);
    for_each(policies::omp, dom, [=](ElemXY const ixy) { view(ixy) += 1; });
    ASSERT_EQ(std::count(storage.begin(), storage.end(), 1), dom.size());
}

static void TestForEachKokkosOneDimension()
{
    DDomX const dom(lbound_x, nelems_x);
    Chunk<int, DDomX, DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ChunkSpan view(storage.span_view());
    for_each(
            policies::parallel_device,
            dom,
            DDC_LAMBDA(ElemX const ix) { view(ix) += 1; });
    int const* const ptr = storage.data();
    int sum;
    Kokkos::parallel_reduce(
            dom.size(),
            KOKKOS_LAMBDA(std::size_t i, int& local_sum) { local_sum += ptr[i]; },
            Kokkos::Sum<int>(sum));
    ASSERT_EQ(sum, dom.size());
}

TEST(ForEachKokkos, OneDimension)
{
    TestForEachKokkosOneDimension();
}

static void TestForEachKokkosTwoDimensions()
{
    DDomXY const dom(lbound_x_y, nelems_x_y);
    Chunk<int, DDomXY, DeviceAllocator<int>> storage(dom);
    Kokkos::deep_copy(storage.allocation_kokkos_view(), 0);
    ChunkSpan view(storage.span_view());
    for_each(
            policies::parallel_device,
            dom,
            DDC_LAMBDA(ElemXY const ixy) { view(ixy) += 1; });
    int const* const ptr = storage.data();
    int sum;
    Kokkos::parallel_reduce(
            dom.size(),
            KOKKOS_LAMBDA(std::size_t i, int& local_sum) { local_sum += ptr[i]; },
            Kokkos::Sum<int>(sum));
    ASSERT_EQ(sum, dom.size());
}

TEST(ForEachKokkos, TwoDimensions)
{
    TestForEachKokkosTwoDimensions();
}
