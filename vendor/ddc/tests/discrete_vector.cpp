// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

struct DDimX;
using DVectX = ddc::DiscreteVector<DDimX>;


struct DDimY;
using DVectY = ddc::DiscreteVector<DDimY>;


struct DDimZ;
using DVectZ = ddc::DiscreteVector<DDimZ>;


using DVectXZ = ddc::DiscreteVector<DDimX, DDimZ>;


using DVectXYZ = ddc::DiscreteVector<DDimX, DDimY, DDimZ>;

} // namespace

TEST(DiscreteVectorTest, ExternalBinaryOperatorPlus)
{
    std::size_t const dv_x = 7;
    std::size_t const dv_y = -2;
    std::size_t const dv_z = 15;
    DVectXYZ dxyz(dv_x, dv_y, dv_z);
    std::size_t const dv2_x = -4;
    std::size_t const dv2_y = 22;
    std::size_t const dv2_z = 8;
    DVectX dx(dv2_x);
    DVectXZ dxz(dv2_x, dv2_z);
    DVectXYZ dxyz2(dv2_x, dv2_y, dv2_z);

    DVectXYZ result1(dxyz + dx);
    DVectXYZ result2(dx + dxyz);
    DVectXYZ result3(dxyz + dxz);
    DVectXYZ result4(dxz + dxyz);
    DVectXYZ result5(dxyz + dxyz2);

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
    DVectXYZ dxyz(dv_x, dv_y, dv_z);
    std::size_t const dv2_x = -4;
    std::size_t const dv2_y = 22;
    std::size_t const dv2_z = 8;
    DVectX dx(dv2_x);
    DVectXZ dxz(dv2_x, dv2_z);
    DVectXYZ dxyz2(dv2_x, dv2_y, dv2_z);

    DVectXYZ result1(dxyz - dx);
    DVectXYZ result2(dx - dxyz);
    DVectXYZ result3(dxyz - dxz);
    DVectXYZ result4(dxz - dxyz);
    DVectXYZ result5(dxyz - dxyz2);

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
