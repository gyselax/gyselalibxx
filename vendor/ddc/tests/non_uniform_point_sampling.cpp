// SPDX-License-Identifier: MIT
#include <array>
#include <sstream>
#include <vector>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#define VALUES_X                                                                                   \
    {                                                                                              \
        ddc::Coordinate<DimX>(0.1), ddc::Coordinate<DimX>(0.2), ddc::Coordinate<DimX>(0.3),        \
                ddc::Coordinate<DimX>(0.4)                                                         \
    }

#define VALUES_Y                                                                                   \
    {                                                                                              \
        ddc::Coordinate<DimY>(0.1), ddc::Coordinate<DimY>(0.2), ddc::Coordinate<DimY>(0.3),        \
                ddc::Coordinate<DimY>(0.4)                                                         \
    }

namespace {

struct DimX;
struct DimY;

using DDimX = ddc::NonUniformPointSampling<DimX>;
using DDimY = ddc::NonUniformPointSampling<DimY>;

std::array<double, 4> const array_points_x VALUES_X;
std::vector<double> const vector_points_x VALUES_X;

std::vector<double> const vector_points_y VALUES_Y;

ddc::DiscreteElement<DDimX> constexpr point_ix(2);
ddc::Coordinate<DimX> constexpr point_rx(0.3);

ddc::DiscreteElement<DDimY> constexpr point_iy(1);
ddc::Coordinate<DimY> constexpr point_ry(0.2);

ddc::DiscreteElement<DDimX, DDimY> constexpr point_ixy(2, 1);
ddc::Coordinate<DimX, DimY> constexpr point_rxy(0.3, 0.2);

} // namespace

TEST(NonUniformPointSamplingTest, ListConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(VALUES_X);
    EXPECT_EQ(ddim_x.size(), 4);
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSamplingTest, ArrayConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(array_points_x);
    EXPECT_EQ(ddim_x.size(), array_points_x.size());
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSamplingTest, VectorConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(vector_points_x);
    EXPECT_EQ(ddim_x.size(), vector_points_x.size());
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSamplingTest, IteratorConstructor)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x(vector_points_x.begin(), vector_points_x.end());
    EXPECT_EQ(ddim_x.size(), vector_points_x.size());
    EXPECT_EQ(ddim_x.coordinate(point_ix), point_rx);
}

TEST(NonUniformPointSampling, Formatting)
{
    DDimX::Impl<Kokkos::HostSpace> ddim_x({ddc::Coordinate<DimX>(0.1), ddc::Coordinate<DimX>(0.4)});
    std::stringstream oss;
    oss << ddim_x;
    EXPECT_EQ(oss.str(), "NonUniformPointSampling(2)");
}

TEST(NonUniformPointSampling, Coordinate)
{
    ddc::init_discrete_space<DDimX>(vector_points_x);
    ddc::init_discrete_space<DDimY>(vector_points_x);
    EXPECT_EQ(coordinate(point_ix), point_rx);
    EXPECT_EQ(coordinate(point_iy), point_ry);
    EXPECT_EQ(coordinate(point_ixy), point_rxy);
}
