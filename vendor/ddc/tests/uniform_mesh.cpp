// SPDX-License-Identifier: MIT
#include <memory>
#include <sstream>

#include <ddc/MCoord>
#include <ddc/RCoord>
#include <ddc/UniformMesh>

#include <gtest/gtest.h>

struct DimX;

using MeshX = UniformMesh<DimX>;

class UniformMeshTest : public ::testing::Test
{
protected:
    RCoord<DimX> origin = RCoord<DimX>(-1.);
    RCoord<DimX> step = RCoord<DimX>(0.5);
    MLength<DimX> npoints = MLength<DimX>(5);
    MCoord<MeshX> lbound = MCoord<MeshX>(0);
    MCoord<MeshX> point_ix = MCoord<MeshX>(2);
    RCoord<DimX> point_rx = RCoord<DimX>(0.);
};

TEST_F(UniformMeshTest, rank)
{
    EXPECT_EQ(MeshX::rank(), 1);
}

TEST_F(UniformMeshTest, constructor)
{
    MeshX mesh_x(origin, step);
    EXPECT_EQ(mesh_x.lbound(), lbound);
    EXPECT_EQ(mesh_x.origin(), origin);
    EXPECT_EQ(mesh_x.step(), step);
    EXPECT_EQ(mesh_x.rmin(), origin);
    EXPECT_EQ(mesh_x.to_real(point_ix), point_rx);
}

TEST_F(UniformMeshTest, constructor2)
{
    MeshX mesh_x(origin, RCoord<DimX>(1.), npoints);
    EXPECT_EQ(mesh_x.lbound(), lbound);
    EXPECT_EQ(mesh_x.origin(), origin);
    EXPECT_EQ(mesh_x.step(), step);
    EXPECT_EQ(mesh_x.rmin(), origin);
    EXPECT_EQ(mesh_x.to_real(point_ix), point_rx);
}

TEST(UniformMesh, formatting)
{
    MeshX mesh_x(RCoord<DimX>(-1.), RCoord<DimX>(0.5));
    std::stringstream oss;
    oss << mesh_x;
    EXPECT_EQ(oss.str(), "UniformMesh( origin=(-1), step=(0.5) )");
}
