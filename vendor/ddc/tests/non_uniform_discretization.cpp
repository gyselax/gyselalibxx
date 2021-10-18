// SPDX-License-Identifier: MIT
#include <array>
#include <sstream>
#include <vector>

#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/NonUniformDiscretization>

#include <gtest/gtest.h>

#define VALUES                                                                                     \
    {                                                                                              \
        Coordinate<DimX>(0.1), Coordinate<DimX>(0.2), Coordinate<DimX>(0.3), Coordinate<DimX>(0.4) \
    }

namespace {

struct DimX;

using DDimX = NonUniformDiscretization<DimX>;

constexpr std::array<double, 4> array_points VALUES;
static std::vector<double> const vector_points VALUES;
constexpr Coordinate<DimX> min_x = Coordinate<DimX>(0.1);
constexpr Coordinate<DimX> max_x = Coordinate<DimX>(0.4);
constexpr DiscreteCoordinate<DDimX> lbound = DiscreteCoordinate<DDimX>(0);
constexpr DiscreteCoordinate<DDimX> ubound = DiscreteCoordinate<DDimX>(3);
constexpr DiscreteCoordinate<DDimX> point_ix = DiscreteCoordinate<DDimX>(2);
constexpr Coordinate<DimX> point_rx = Coordinate<DimX>(0.3);

} // namespace

TEST(NonUniformDiscretizationTest, rank)
{
    EXPECT_EQ(DDimX::rank(), 1);
}

TEST(NonUniformDiscretizationTest, list_constructor)
{
    DDimX ddim_x(VALUES);
    EXPECT_EQ(ddim_x.size(), 4);
    EXPECT_EQ(ddim_x.lbound(), lbound);
    EXPECT_EQ(ddim_x.ubound(), ubound);
    EXPECT_EQ(ddim_x.rmin(), min_x);
    EXPECT_EQ(ddim_x.rmax(), max_x);
    EXPECT_EQ(ddim_x.rlength(), max_x - min_x);
    EXPECT_EQ(ddim_x.to_real(point_ix), point_rx);
}

TEST(NonUniformDiscretizationTest, array_constructor)
{
    DDimX ddim_x(array_points);
    EXPECT_EQ(ddim_x.size(), array_points.size());
    EXPECT_EQ(ddim_x.lbound(), lbound);
    EXPECT_EQ(ddim_x.ubound(), ubound);
    EXPECT_EQ(ddim_x.rmin(), min_x);
    EXPECT_EQ(ddim_x.rmax(), max_x);
    EXPECT_EQ(ddim_x.rlength(), max_x - min_x);
    EXPECT_EQ(ddim_x.to_real(point_ix), point_rx);
}

TEST(NonUniformDiscretizationTest, vector_constructor)
{
    DDimX ddim_x(vector_points);
    EXPECT_EQ(ddim_x.size(), vector_points.size());
    EXPECT_EQ(ddim_x.lbound(), lbound);
    EXPECT_EQ(ddim_x.ubound(), ubound);
    EXPECT_EQ(ddim_x.rmin(), min_x);
    EXPECT_EQ(ddim_x.rmax(), max_x);
    EXPECT_EQ(ddim_x.rlength(), max_x - min_x);
    EXPECT_EQ(ddim_x.to_real(point_ix), point_rx);
}

TEST(NonUniformDiscretizationTest, iterator_constructor)
{
    DDimX ddim_x(vector_points.begin(), vector_points.end());
    EXPECT_EQ(ddim_x.size(), vector_points.size());
    EXPECT_EQ(ddim_x.lbound(), lbound);
    EXPECT_EQ(ddim_x.ubound(), ubound);
    EXPECT_EQ(ddim_x.rmin(), min_x);
    EXPECT_EQ(ddim_x.rmax(), max_x);
    EXPECT_EQ(ddim_x.rlength(), max_x - min_x);
    EXPECT_EQ(ddim_x.to_real(point_ix), point_rx);
}

TEST(NonUniformDiscretization, formatting)
{
    DDimX ddim_x({Coordinate<DimX>(0.1), Coordinate<DimX>(0.4)});
    std::stringstream oss;
    oss << ddim_x;
    EXPECT_EQ(oss.str(), "NonUniformDiscretization( (0.1), ..., (0.4) )");
}
