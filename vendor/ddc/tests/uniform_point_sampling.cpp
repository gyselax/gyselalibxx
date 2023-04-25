// SPDX-License-Identifier: MIT
#include <memory>
#include <sstream>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

struct DimX;

using DDimX = ddc::UniformPointSampling<DimX>;

class UniformPointSamplingTest : public ::testing::Test
{
protected:
    ddc::Coordinate<DimX> origin = ddc::Coordinate<DimX>(-1.);
    double step = 0.5;
    ddc::DiscreteVector<DDimX> npoints = ddc::DiscreteVector<DDimX>(5);
    ddc::DiscreteElement<DDimX> lbound = ddc::DiscreteElement<DDimX>(0);
    ddc::DiscreteElement<DDimX> point_ix = ddc::DiscreteElement<DDimX>(2);
    ddc::Coordinate<DimX> point_rx = ddc::Coordinate<DimX>(0.);
};

TEST_F(UniformPointSamplingTest, Constructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(origin, step);
    EXPECT_EQ(ddim_x.origin(), origin);
    EXPECT_EQ(ddim_x.step(), step);
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(UniformPointSampling, Formatting)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(ddc::Coordinate<DimX>(-1.), 0.5);
    std::stringstream oss;
    oss << ddim_x;
    EXPECT_EQ(oss.str(), "UniformPointSampling( origin=(-1), step=0.5 )");
}
