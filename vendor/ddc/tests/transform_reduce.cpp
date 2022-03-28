// SPDX-License-Identifier: MIT

#include <ddc/Chunk>
#include <ddc/reducer.hpp>
#include <ddc/transform_reduce>

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

TEST(TransformReduceSerial, OneDimension)
{
    DDomX const dom(lbound_x, nelems_x);
    Chunk<int, DDomX> chunk(dom);
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
    Chunk<int, DDomXY> chunk(dom);
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
    Chunk<int, DDomX> chunk(dom);
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
    Chunk<int, DDomXY> chunk(dom);
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
