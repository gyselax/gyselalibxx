#include <gtest/gtest.h>

#include "mesh_product.h"
#include "nonuniformmesh.h"
#include "uniformmesh.h"

class DimX;
class DimVx;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<DimX>;
using RCoordVx = RCoord<DimVx>;
using MCoordVx = MCoord<DimVx>;
using MeshX = UniformMesh<DimX>;
using MeshVx = NonUniformMesh<DimVx>;

using MeshXVx = MeshProduct<MeshX, MeshVx>;
using MCoordXVx = MeshXVx::mcoord_type;
using RCoordXVx = MeshXVx::rcoord_type;

TEST(MeshProduct, constructor)
{
    MeshX mesh_x(2., 0.1);

    std::array points_vx {-1., 0., 2., 4.};
    MeshVx mesh_vx(points_vx, points_vx.size());

    MeshProduct mesh_x_vx(mesh_x, mesh_vx);
    EXPECT_EQ(MeshXVx::rank(), MeshX::rank() + MeshVx::rank());
}

TEST(MeshProduct, submesh)
{
    MeshX mesh_x(2., 0.1);

    std::array points_vx {-1., 0., 2., 4.};
    MeshVx mesh_vx(points_vx, points_vx.size());

    MeshProduct mesh_x_vx(mesh_x, mesh_vx);
    auto&& submesh = mesh_x_vx.submesh(2, std::experimental::all);
    EXPECT_EQ(1, submesh.rank());
    EXPECT_EQ(RCoordXVx(2.2, 0.), submesh.to_real(MCoordXVx(0, 1)));
}

TEST(MeshProduct, to_real)
{
    MeshX mesh_x(0., 0.1);

    std::array points_vx {-1., 0., 2., 4.};
    MeshVx mesh_vx(points_vx, points_vx.size());

    MeshProduct mesh_x_vx(mesh_x, mesh_vx);

    for (std::size_t ix = 0; ix < 5; ++ix) {
        for (std::size_t ivx = 0; ivx < points_vx.size(); ++ivx) {
            EXPECT_EQ(
                    RCoordXVx(mesh_x.to_real(ix), points_vx[ivx]),
                    mesh_x_vx.to_real(MCoordXVx(ix, ivx)));
        }
    }
}
