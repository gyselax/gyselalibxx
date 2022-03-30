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
static DDomXY constexpr dom_x_y(lbound_x_y, nelems_x_y);

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
