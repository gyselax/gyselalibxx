// SPDX-License-Identifier: MIT
#include <memory>
#include <utility>

#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
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

using MCoordXVx = MCoord<MeshX, MeshVx>;
using RCoordXVx = RCoord<DimX, DimVx>;

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
