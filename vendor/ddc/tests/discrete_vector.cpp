// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <cstddef>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace DDC_HIP_5_7_ANONYMOUS_NAMESPACE_WORKAROUND(DISCRETE_VECTOR_CPP)
{
    struct DDimX
    {
    };
    using DVectX = ddc::DiscreteVector<DDimX>;


    struct DDimY
    {
    };
    using DVectY = ddc::DiscreteVector<DDimY>;


    struct DDimZ
    {
    };
    using DVectZ = ddc::DiscreteVector<DDimZ>;


    using DVectXZ = ddc::DiscreteVector<DDimX, DDimZ>;


    using DVectXYZ = ddc::DiscreteVector<DDimX, DDimY, DDimZ>;

} // namespace )

TEST(DiscreteVectorXYZTest, ConstructorFromDiscreteVectors)
{
    std::size_t const dv_x = 7;
    std::size_t const dv_y = 13;
    std::size_t const dv_z = 4;
    DVectXZ const ixz(dv_x, dv_z);
    DVectY const iy(dv_y);
    DVectXYZ const ixyz(ixz, iy);
    EXPECT_EQ(ixyz.get<DDimX>(), dv_x);
    EXPECT_EQ(ixyz.get<DDimY>(), dv_y);
    EXPECT_EQ(ixyz.get<DDimZ>(), dv_z);
}

TEST(DiscreteVectorXTest, PreIncrement)
{
    DVectX const ix0(3);
    DVectX ix1(ix0);
    DVectX const ix2 = ++ix1;
    EXPECT_EQ(ix1, ix0 + 1);
    EXPECT_EQ(ix2, ix0 + 1);
}

TEST(DiscreteVectorXTest, PostIncrement)
{
    DVectX const ix0(3);
    DVectX ix1(ix0);
    DVectX const ix2 = ix1++;
    EXPECT_EQ(ix1, ix0 + 1);
    EXPECT_EQ(ix2, ix0);
}

TEST(DiscreteVectorXTest, PreDecrement)
{
    DVectX const ix0(3);
    DVectX ix1(ix0);
    DVectX const ix2 = --ix1;
    EXPECT_EQ(ix1, ix0 - 1);
    EXPECT_EQ(ix2, ix0 - 1);
}

TEST(DiscreteVectorXTest, PostDecrement)
{
    DVectX const ix0(3);
    DVectX ix1(ix0);
    DVectX const ix2 = ix1--;
    EXPECT_EQ(ix1, ix0 - 1);
    EXPECT_EQ(ix2, ix0);
}

TEST(DiscreteVectorTest, ExternalBinaryOperatorPlus)
{
    std::size_t const dv_x = 7;
    std::size_t const dv_y = -2;
    std::size_t const dv_z = 15;
    DVectXYZ const dxyz(dv_x, dv_y, dv_z);
    std::size_t const dv2_x = -4;
    std::size_t const dv2_y = 22;
    std::size_t const dv2_z = 8;
    DVectX const dx(dv2_x);
    DVectXZ const dxz(dv2_x, dv2_z);
    DVectXYZ const dxyz2(dv2_x, dv2_y, dv2_z);

    DVectXYZ const result1(dxyz + dx);
    DVectXYZ const result2(dx + dxyz);
    DVectXYZ const result3(dxyz + dxz);
    DVectXYZ const result4(dxz + dxyz);
    DVectXYZ const result5(dxyz + dxyz2);

    EXPECT_EQ(ddc::get<DDimX>(result1), dv_x + dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result1), dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result1), dv_z);

    EXPECT_EQ(ddc::get<DDimX>(result2), dv_x + dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result2), dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result2), dv_z);

    EXPECT_EQ(ddc::get<DDimX>(result3), dv_x + dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result3), dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result3), dv_z + dv2_z);

    EXPECT_EQ(ddc::get<DDimX>(result4), dv_x + dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result4), dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result4), dv_z + dv2_z);

    EXPECT_EQ(ddc::get<DDimX>(result5), dv_x + dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result5), dv_y + dv2_y);
    EXPECT_EQ(ddc::get<DDimZ>(result5), dv_z + dv2_z);
}

TEST(DiscreteVectorTest, ExternalBinaryOperatorMinus)
{
    std::size_t const dv_x = 7;
    std::size_t const dv_y = -2;
    std::size_t const dv_z = 15;
    DVectXYZ const dxyz(dv_x, dv_y, dv_z);
    std::size_t const dv2_x = -4;
    std::size_t const dv2_y = 22;
    std::size_t const dv2_z = 8;
    DVectX const dx(dv2_x);
    DVectXZ const dxz(dv2_x, dv2_z);
    DVectXYZ const dxyz2(dv2_x, dv2_y, dv2_z);

    DVectXYZ const result1(dxyz - dx);
    DVectXYZ const result2(dx - dxyz);
    DVectXYZ const result3(dxyz - dxz);
    DVectXYZ const result4(dxz - dxyz);
    DVectXYZ const result5(dxyz - dxyz2);

    EXPECT_EQ(ddc::get<DDimX>(result1), dv_x - dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result1), dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result1), dv_z);

    EXPECT_EQ(ddc::get<DDimX>(result2), -dv_x + dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result2), -dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result2), -dv_z);

    EXPECT_EQ(ddc::get<DDimX>(result3), dv_x - dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result3), dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result3), dv_z - dv2_z);

    EXPECT_EQ(ddc::get<DDimX>(result4), -dv_x + dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result4), -dv_y);
    EXPECT_EQ(ddc::get<DDimZ>(result4), -dv_z + dv2_z);

    EXPECT_EQ(ddc::get<DDimX>(result5), dv_x - dv2_x);
    EXPECT_EQ(ddc::get<DDimY>(result5), dv_y - dv2_y);
    EXPECT_EQ(ddc::get<DDimZ>(result5), dv_z - dv2_z);
}
