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
using MeshX = NonUniformMesh<DimX>;
using MeshVx = NonUniformMesh<DimVx>;

using MeshXVx = MeshProduct<MeshX, MeshVx>;
using MCoordXVx = MCoord<DimX, DimVx>;
using RCoordXVx = RCoord<DimX, DimVx>;

TEST(MeshProduct, constructor)
{
    std::array points_x {0., 1., 2.};
    MeshX mesh_x(points_x, points_x.size());

    std::array points_vx {-1., 0., 2., 4.};
    MeshVx mesh_vx(points_vx, points_vx.size());

    MeshProduct mesh_x_vx(mesh_x, mesh_vx);
    EXPECT_EQ(MeshXVx::rank(), MeshX::rank() + MeshVx::rank());
}

TEST(MeshProduct, to_real)
{
    std::array points_x {0., 1., 2.};
    MeshX mesh_x(points_x, points_x.size());

    std::array points_vx {-1., 0., 2., 4.};
    MeshVx mesh_vx(points_vx, points_vx.size());

    MeshProduct mesh_x_vx(mesh_x, mesh_vx);

    for (std::size_t ix = 0; ix < points_x.size(); ++ix) {
        for (std::size_t ivx = 0; ivx < points_vx.size(); ++ivx) {
            EXPECT_EQ(
                    RCoordXVx(points_x[ix], points_vx[ivx]),
                    mesh_x_vx.to_real(MCoordXVx(ix, ivx)));
        }
    }
}
