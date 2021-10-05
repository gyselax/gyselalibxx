// SPDX-License-Identifier: MIT
#include <array>
#include <iosfwd>
#include <memory>

#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>
#include <ddc/detail/product_mesh.hpp>

#include <gtest/gtest.h>

class DimX;
class DimVx;

using RCoordX = RCoord<DimX>;
using RCoordVx = RCoord<DimVx>;
using MeshX = UniformMesh<DimX>;
using MeshVx = NonUniformMesh<DimVx>;

using MeshXVx = detail::ProductMesh<MeshX, MeshVx>;
using MCoordXVx = MCoord<MeshX, MeshVx>;
using RCoordXVx = RCoord<DimX, DimVx>;

class ProductMeshTest : public ::testing::Test
{
protected:
    MeshX mesh_x = MeshX(2., 0.1);
    std::array<RCoordVx, 4> points_vx {-1., 0., 2., 4.};
    MeshVx mesh_vx = MeshVx(points_vx);
    MeshXVx mesh_x_vx {mesh_x, mesh_vx};
};

TEST_F(ProductMeshTest, constructor)
{
    EXPECT_EQ(MeshXVx::rank(), MeshX::rank() + MeshVx::rank());
    EXPECT_EQ(mesh_x_vx.to_real(MCoordXVx(0, 0)), RCoordXVx(2., -1.));
}

TEST_F(ProductMeshTest, accessor)
{
    EXPECT_EQ(mesh_x_vx.get<MeshX>(), mesh_x);
}

TEST_F(ProductMeshTest, submesh)
{
    auto&& selection = select<MeshVx>(mesh_x_vx);
    EXPECT_EQ(1, selection.rank());
    EXPECT_EQ(RCoordVx(0.), selection.to_real(MCoord<MeshVx>(1)));
}

TEST_F(ProductMeshTest, conversion)
{
    constexpr static MeshX mesh_x(2., 0.1);
    constexpr detail::ProductMesh product_mesh_x(mesh_x);
    MeshX const& mesh_x_ref = get<MeshX>(product_mesh_x);
    double step = mesh_x_ref.step();
    EXPECT_EQ(0.1, step);
}

TEST_F(ProductMeshTest, to_real)
{
    for (std::size_t ix = 0; ix < 5; ++ix) {
        for (std::size_t ivx = 0; ivx < points_vx.size(); ++ivx) {
            EXPECT_EQ(
                    RCoordXVx(mesh_x.to_real(ix), points_vx[ivx]),
                    mesh_x_vx.to_real(MCoordXVx(ix, ivx)));
        }
    }
}
