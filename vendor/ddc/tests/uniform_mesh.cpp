#include <memory>

#include <gtest/gtest.h>

#include "gtest/gtest_pred_impl.h"

#include "mcoord.h"
#include "rcoord.h"
#include "uniform_mesh.h"

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
