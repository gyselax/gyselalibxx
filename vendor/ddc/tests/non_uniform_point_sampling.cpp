// SPDX-License-Identifier: MIT
#include <array>
#include <sstream>
#include <vector>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#define VALUES                                                                                     \
    {                                                                                              \
        Coordinate<DimX>(0.1), Coordinate<DimX>(0.2), Coordinate<DimX>(0.3), Coordinate<DimX>(0.4) \
    }

namespace {
struct DimX;

using DDimX = NonUniformPointSampling<DimX>;

static std::array<double, 4> const array_points VALUES;
static std::vector<double> const vector_points VALUES;
DiscreteElement<DDimX> constexpr point_ix = DiscreteElement<DDimX>(2);
Coordinate<DimX> constexpr point_rx = Coordinate<DimX>(0.3);

} // namespace

TEST(NonUniformPointSamplingTest, ListConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(VALUES);
    EXPECT_EQ(ddim_x.size(), 4);
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSamplingTest, ArrayConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(array_points);
    EXPECT_EQ(ddim_x.size(), array_points.size());
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSamplingTest, VectorConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(vector_points);
    EXPECT_EQ(ddim_x.size(), vector_points.size());
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSamplingTest, IteratorConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(vector_points.begin(), vector_points.end());
    EXPECT_EQ(ddim_x.size(), vector_points.size());
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSampling, Formatting)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x({Coordinate<DimX>(0.1), Coordinate<DimX>(0.4)});
    std::stringstream oss;
    oss << ddim_x;
    EXPECT_EQ(oss.str(), "NonUniformPointSampling(2)");
}
