#include <gtest/gtest.h>

#include "non_uniform_mesh.h"
#include "product_mesh.h"
#include "single_mesh.h"
#include "uniform_mesh.h"

class DimX;
class DimVx;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<DimX>;
using RCoordVx = RCoord<DimVx>;
using MCoordVx = MCoord<DimVx>;
using MeshX = SingleMesh<DimX>;
using MeshVx = NonUniformMesh<DimVx>;

using MeshXVx = ProductMesh<MeshX, MeshVx>;
using MCoordXVx = MeshXVx::mcoord_type;
using RCoordXVx = MeshXVx::rcoord_type;

TEST(SingleMesh, class_size)
{
    EXPECT_EQ(sizeof(MeshX), sizeof(double));
}

TEST(SingleMesh, constructor)
{
    constexpr RCoordX x = 1.;

    SingleMesh<DimX> mesh_x(x);

    EXPECT_EQ(mesh_x.to_real(MCoordX(0)), x);
}

TEST(SingleMesh, product)
{
    SingleMesh<DimX> mesh_x(RCoordX(1.));

    std::array points_vx {-1., 0., 2., 4.};
    MeshVx mesh_vx(points_vx, points_vx.size());

    ProductMesh mesh_x_vx(mesh_x, mesh_vx);
    EXPECT_EQ(MeshXVx::rank(), MeshX::rank() + MeshVx::rank());
}
