// SPDX-License-Identifier: MIT
#include <array>
#include <sstream>
#include <vector>

#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

#include <gtest/gtest.h>

struct DimX;

using MeshX = NonUniformMesh<DimX>;

#define VALUES                                                                                     \
    {                                                                                              \
        0.1, 0.2, 0.3, 0.4                                                                         \
    }

class NonUniformMeshTest : public testing::Test
{
protected:
    std::array<double, 4> array_points VALUES;
    std::vector<double> vector_points VALUES;
    RCoord<DimX> min_x = 0.1;
    RCoord<DimX> max_x = 0.4;
    MCoord<MeshX> lbound = 0;
    MCoord<MeshX> ubound = 3;
    MCoord<MeshX> point_ix = 2;
    RCoord<DimX> point_rx = 0.3;
};

TEST_F(NonUniformMeshTest, rank)
{
    EXPECT_EQ(MeshX::rank(), 1);
}

TEST_F(NonUniformMeshTest, list_constructor)
{
    MeshX mesh_x(VALUES);
    EXPECT_EQ(mesh_x.size(), 4);
    EXPECT_EQ(mesh_x.lbound(), lbound);
    EXPECT_EQ(mesh_x.ubound(), ubound);
    EXPECT_EQ(mesh_x.rmin(), min_x);
    EXPECT_EQ(mesh_x.rmax(), max_x);
    EXPECT_EQ(mesh_x.rlength(), max_x - min_x);
    EXPECT_EQ(mesh_x.to_real(point_ix), point_rx);
}

TEST_F(NonUniformMeshTest, array_constructor)
{
    MeshX mesh_x(array_points);
    EXPECT_EQ(mesh_x.size(), array_points.size());
    EXPECT_EQ(mesh_x.lbound(), lbound);
    EXPECT_EQ(mesh_x.ubound(), ubound);
    EXPECT_EQ(mesh_x.rmin(), min_x);
    EXPECT_EQ(mesh_x.rmax(), max_x);
    EXPECT_EQ(mesh_x.rlength(), max_x - min_x);
    EXPECT_EQ(mesh_x.to_real(point_ix), point_rx);
}

TEST_F(NonUniformMeshTest, vector_constructor)
{
    MeshX mesh_x(vector_points);
    EXPECT_EQ(mesh_x.size(), vector_points.size());
    EXPECT_EQ(mesh_x.lbound(), lbound);
    EXPECT_EQ(mesh_x.ubound(), ubound);
    EXPECT_EQ(mesh_x.rmin(), min_x);
    EXPECT_EQ(mesh_x.rmax(), max_x);
    EXPECT_EQ(mesh_x.rlength(), max_x - min_x);
    EXPECT_EQ(mesh_x.to_real(point_ix), point_rx);
}

TEST_F(NonUniformMeshTest, iterator_constructor)
{
    MeshX mesh_x(vector_points.begin(), vector_points.end());
    EXPECT_EQ(mesh_x.size(), vector_points.size());
    EXPECT_EQ(mesh_x.lbound(), lbound);
    EXPECT_EQ(mesh_x.ubound(), ubound);
    EXPECT_EQ(mesh_x.rmin(), min_x);
    EXPECT_EQ(mesh_x.rmax(), max_x);
    EXPECT_EQ(mesh_x.rlength(), max_x - min_x);
    EXPECT_EQ(mesh_x.to_real(point_ix), point_rx);
}

TEST(NonUniformMesh, formatting)
{
    MeshX mesh_x({0.1, 0.4});
    std::stringstream oss;
    oss << mesh_x;
    EXPECT_EQ(oss.str(), "NonUniformMesh( (0.1), ..., (0.4) )");
}
