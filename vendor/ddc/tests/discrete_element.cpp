// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

struct DDimX;
using ElemX = DiscreteCoordinate<DDimX>;
using DVectX = DiscreteVector<DDimX>;


struct DDimY;
using ElemY = DiscreteCoordinate<DDimY>;
using DVectY = DiscreteVector<DDimY>;


using ElemXY = DiscreteCoordinate<DDimX, DDimY>;
using DVectXY = DiscreteVector<DDimX, DDimY>;


using ElemYX = DiscreteCoordinate<DDimY, DDimX>;
using DVectYX = DiscreteVector<DDimY, DDimX>;

} // namespace


TEST(DiscreteElementXTest, RightExternalBinaryOperatorPlus)
{
    std::size_t const uid_x = 7;
    ElemX const ix(uid_x);
    std::ptrdiff_t const dv_x = -2;
    ElemX const ix2 = ix + DVectX(dv_x);
    ASSERT_EQ(ix2.uid<DDimX>(), uid_x + dv_x);
}

TEST(DiscreteElementXTest, RightExternalBinaryOperatorMinus)
{
    std::size_t const uid_x = 7;
    ElemX const ix(uid_x);
    std::ptrdiff_t const dv_x = -2;
    ElemX const ixy2 = ix - dv_x;
    ASSERT_EQ(ixy2.uid<DDimX>(), uid_x - dv_x);
}

TEST(DiscreteElementXTest, BinaryOperatorMinus)
{
    std::size_t const uid_x = 7;
    ElemX const ix(uid_x);
    std::ptrdiff_t const dv_x = -2;
    ElemX const ix2 = ix + dv_x;
    DVectX dv2_x = ix2 - ix;
    ASSERT_EQ(get<DDimX>(dv2_x), dv_x);
}

TEST(DiscreteElementXYTest, ValueConstructor)
{
    ElemXY const ixy {};
    ASSERT_EQ(ixy.uid<DDimX>(), std::size_t());
    ASSERT_EQ(ixy.uid<DDimY>(), std::size_t());
}

TEST(DiscreteElementXYTest, UntaggedConstructor)
{
    std::size_t const uid_x = 7;
    std::size_t const uid_y = 13;
    ElemXY const ixy(uid_x, uid_y);
    ASSERT_EQ(ixy.uid<DDimX>(), uid_x);
    ASSERT_EQ(ixy.uid<DDimY>(), uid_y);
}

TEST(DiscreteElementXYTest, RightExternalBinaryOperatorPlus)
{
    std::size_t const uid_x = 7;
    std::size_t const uid_y = 13;
    ElemXY const ixy(uid_x, uid_y);
    std::ptrdiff_t const dv_x = -2;
    std::ptrdiff_t const dv_y = +3;
    ElemXY const ixy2 = ixy + DVectXY(dv_x, dv_y);
    ASSERT_EQ(ixy2.uid<DDimX>(), uid_x + dv_x);
    ASSERT_EQ(ixy2.uid<DDimY>(), uid_y + dv_y);
}

TEST(DiscreteElementXYTest, RightExternalBinaryOperatorMinus)
{
    std::size_t const uid_x = 7;
    std::size_t const uid_y = 13;
    ElemXY const ixy(uid_x, uid_y);
    std::ptrdiff_t const dv_x = -2;
    std::ptrdiff_t const dv_y = +3;
    ElemXY const ixy2 = ixy - DVectXY(dv_x, dv_y);
    ASSERT_EQ(ixy2.uid<DDimX>(), uid_x - dv_x);
    ASSERT_EQ(ixy2.uid<DDimY>(), uid_y - dv_y);
}

TEST(DiscreteElementXYTest, BinaryOperatorMinus)
{
    std::size_t const uid_x = 7;
    std::size_t const uid_y = 13;
    ElemXY const ixy(uid_x, uid_y);
    std::ptrdiff_t const dv_x = -2;
    std::ptrdiff_t const dv_y = +3;
    ElemXY const ixy2 = ixy + DVectXY(dv_x, dv_y);
    DVectXY dv_xy = ixy2 - ixy;
    ASSERT_EQ(get<DDimX>(dv_xy), dv_x);
    ASSERT_EQ(get<DDimY>(dv_xy), dv_y);
}
