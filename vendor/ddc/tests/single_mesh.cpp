#include <memory>
#include <utility>

#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
#include <ddc/ProductMesh>
#include <ddc/RCoord>
#include <ddc/SingleMesh>

#include <gtest/gtest.h>

class DimX;
class DimVx;

using RCoordX = RCoord<DimX>;
using RCoordVx = RCoord<DimVx>;
using MeshX = SingleMesh<DimX>;
using MeshVx = NonUniformMesh<DimVx>;

using MCoordX = MCoord<MeshX>;
using MCoordVx = MCoord<MeshVx>;

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
    MeshVx mesh_vx(points_vx);

    ProductMesh mesh_x_vx(mesh_x, mesh_vx);
    EXPECT_EQ(MeshXVx::rank(), MeshX::rank() + MeshVx::rank());
}
