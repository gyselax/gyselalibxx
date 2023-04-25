// SPDX-License-Identifier: MIT
#include <memory>
#include <sstream>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace {

struct DimX;

using DDimX = ddc::UniformPointSampling<DimX>;

static ddc::Coordinate<DimX> constexpr origin(-1.);
static double constexpr step = 0.5;
static ddc::DiscreteElement<DDimX> constexpr point_ix(2);
static ddc::Coordinate<DimX> constexpr point_rx(0.);

} // namespace

TEST(UniformPointSamplingTest, Constructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(origin, step);
    EXPECT_EQ(ddim_x.origin(), origin);
    EXPECT_EQ(ddim_x.step(), step);
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(UniformPointSampling, Formatting)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(origin, step);
    std::stringstream oss;
    oss << ddim_x;
    EXPECT_EQ(oss.str(), "UniformPointSampling( origin=(-1), step=0.5 )");
}
