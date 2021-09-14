#include <array>
#include <vector>

#include <ddc/mcoord.h>
#include <ddc/non_uniform_mesh.h>
#include <ddc/rcoord.h>
#include <ddc/taggedvector.h>

#include <gtest/gtest.h>

struct DimX;

using MeshX = NonUniformMesh<DimX>;

class NonUniformMeshTest : public testing::Test
{
protected:
    std::array<double, 4> array_points {0.1, 0.2, 0.3, 0.4};
    std::vector<double> vector_points {0.1, 0.2, 0.3, 0.4};
};

TEST_F(NonUniformMeshTest, list_constructor)
{
    MeshX mesh_x({0.1, 0.2, 0.3, 0.4});
    EXPECT_EQ(mesh_x.rmin(), 0.1);
    EXPECT_EQ(mesh_x.rmax(), 0.4);
    EXPECT_EQ(mesh_x.to_real(MCoord<MeshX>(2)), RCoord<DimX>(0.3));
}

TEST_F(NonUniformMeshTest, array_constructor)
{
    MeshX mesh_x(array_points);
    EXPECT_EQ(mesh_x.rmin(), 0.1);
    EXPECT_EQ(mesh_x.rmax(), 0.4);
    EXPECT_EQ(mesh_x.to_real(MCoord<MeshX>(2)), RCoord<DimX>(0.3));
}

TEST_F(NonUniformMeshTest, vector_constructor)
{
    MeshX mesh_x(vector_points);
    EXPECT_EQ(mesh_x.rmin(), 0.1);
    EXPECT_EQ(mesh_x.rmax(), 0.4);
    EXPECT_EQ(mesh_x.to_real(MCoord<MeshX>(2)), RCoord<DimX>(0.3));
}
