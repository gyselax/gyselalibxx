// SPDX-License-Identifier: MIT
#include <memory>

#include <ddc/MCoord>
#include <ddc/RCoord>
#include <ddc/UniformMesh>

#include <gtest/gtest.h>

struct DimX;

using MeshX = UniformMesh<DimX>;

TEST(UniformMeshNewTest, constructor)
{
    MeshX mesh_x(-1., 0.5);
    EXPECT_EQ(mesh_x.origin(), RCoord<DimX>(-1.));
    EXPECT_EQ(mesh_x.step(), RCoord<DimX>(0.5));
    EXPECT_EQ(mesh_x.to_real(MCoord<MeshX>(2)), RCoord<DimX>(0.));
}

TEST(UniformMeshNewTest, constructor2)
{
    MeshX mesh_x(-1., 1., 5);
    EXPECT_EQ(mesh_x.origin(), RCoord<DimX>(-1.));
    EXPECT_EQ(mesh_x.step(), RCoord<DimX>(0.5));
    EXPECT_EQ(mesh_x.to_real(MCoord<MeshX>(2)), RCoord<DimX>(0.));
}
