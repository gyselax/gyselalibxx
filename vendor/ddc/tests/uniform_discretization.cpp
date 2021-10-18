// SPDX-License-Identifier: MIT
#include <memory>
#include <sstream>

#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/UniformDiscretization>

#include <gtest/gtest.h>

struct DimX;

using DDimX = UniformDiscretization<DimX>;

class UniformDiscretizationTest : public ::testing::Test
{
protected:
    Coordinate<DimX> origin = Coordinate<DimX>(-1.);
    Coordinate<DimX> step = Coordinate<DimX>(0.5);
    DiscreteVector<DimX> npoints = DiscreteVector<DimX>(5);
    DiscreteCoordinate<DDimX> lbound = DiscreteCoordinate<DDimX>(0);
    DiscreteCoordinate<DDimX> point_ix = DiscreteCoordinate<DDimX>(2);
    Coordinate<DimX> point_rx = Coordinate<DimX>(0.);
};

TEST_F(UniformDiscretizationTest, rank)
{
    EXPECT_EQ(DDimX::rank(), 1);
}

TEST_F(UniformDiscretizationTest, constructor)
{
    DDimX ddim_x(origin, step);
    EXPECT_EQ(ddim_x.lbound(), lbound);
    EXPECT_EQ(ddim_x.origin(), origin);
    EXPECT_EQ(ddim_x.step(), step);
    EXPECT_EQ(ddim_x.rmin(), origin);
    EXPECT_EQ(ddim_x.to_real(point_ix), point_rx);
}

TEST_F(UniformDiscretizationTest, constructor2)
{
    DDimX ddim_x(origin, Coordinate<DimX>(1.), npoints);
    EXPECT_EQ(ddim_x.lbound(), lbound);
    EXPECT_EQ(ddim_x.origin(), origin);
    EXPECT_EQ(ddim_x.step(), step);
    EXPECT_EQ(ddim_x.rmin(), origin);
    EXPECT_EQ(ddim_x.to_real(point_ix), point_rx);
}

TEST(UniformDiscretization, formatting)
{
    DDimX ddim_x(Coordinate<DimX>(-1.), Coordinate<DimX>(0.5));
    std::stringstream oss;
    oss << ddim_x;
    EXPECT_EQ(oss.str(), "UniformDiscretization( origin=(-1), step=(0.5) )");
}
