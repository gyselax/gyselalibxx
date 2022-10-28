// SPDX-License-Identifier: MIT
#include <memory>
#include <sstream>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

struct DimX;

using DDimX = UniformPointSampling<DimX>;

class UniformPointSamplingTest : public ::testing::Test
{
protected:
    Coordinate<DimX> origin = Coordinate<DimX>(-1.);
    Coordinate<DimX> step = Coordinate<DimX>(0.5);
    DiscreteVector<DDimX> npoints = DiscreteVector<DDimX>(5);
    DiscreteElement<DDimX> lbound = DiscreteElement<DDimX>(0);
    DiscreteElement<DDimX> point_ix = DiscreteElement<DDimX>(2);
    Coordinate<DimX> point_rx = Coordinate<DimX>(0.);
};

TEST_F(UniformPointSamplingTest, Constructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(origin, step);
    EXPECT_EQ(ddim_x.origin(), origin);
    EXPECT_EQ(ddim_x.step(), step);
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST_F(UniformPointSamplingTest, Constructor2)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(origin, Coordinate<DimX>(1.), npoints);
    EXPECT_EQ(ddim_x.origin(), origin);
    EXPECT_EQ(ddim_x.step(), step);
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(UniformPointSampling, Formatting)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(Coordinate<DimX>(-1.), Coordinate<DimX>(0.5));
    std::stringstream oss;
    oss << ddim_x;
    EXPECT_EQ(oss.str(), "UniformPointSampling( origin=(-1), step=(0.5) )");
}
