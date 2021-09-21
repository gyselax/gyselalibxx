// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
#include <ddc/ProductMDomain>
#include <ddc/ProductMesh>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>

#include <gtest/gtest.h>

class DimX;
class DimVx;

using MeshX = UniformMesh<DimX>;
using MeshVx = NonUniformMesh<DimVx>;

class ProductMDomainTest : public ::testing::Test
{
protected:
    std::size_t npoints_x = 11;
    std::size_t npoints_vx = 4;
    MeshX mesh_x = MeshX(2., 3., npoints_x);
    MeshVx mesh_vx = MeshVx({-1., 0., 2., 4.});
    ProductMesh<MeshX, MeshVx> mesh_x_vx = ProductMesh(mesh_x, mesh_vx);
    ProductMDomain<MeshX, MeshVx> domain_x_vx = ProductMDomain(
            mesh_x_vx,
            MCoord<MeshX, MeshVx>(0, 0),
            MLength<MeshX, MeshVx>(npoints_x, npoints_vx));
};

TEST_F(ProductMDomainTest, constructor)
{
    EXPECT_EQ(domain_x_vx.extents(), (MCoord<MeshX, MeshVx>(npoints_x, npoints_vx)));
    EXPECT_EQ(domain_x_vx.front(), (MCoord<MeshX, MeshVx>(0, 0)));
    EXPECT_EQ(domain_x_vx.back(), (MCoord<MeshX, MeshVx>(npoints_x - 1, npoints_vx - 1)));
    EXPECT_EQ(domain_x_vx.size(), npoints_x * npoints_vx);
}

TEST_F(ProductMDomainTest, rmin_rmax)
{
    EXPECT_EQ(domain_x_vx.rmin(), (RCoord<DimX, DimVx>(2., -1.)));
    EXPECT_EQ(domain_x_vx.rmax(), (RCoord<DimX, DimVx>(3., 4.)));
}

TEST_F(ProductMDomainTest, subdomain)
{
    ProductMDomain subdomain_x(ProductMesh<MeshX>(mesh_x), MCoord<MeshX>(1), MLength<MeshX>(1));
    auto subdomain = domain_x_vx.restrict(subdomain_x);
    EXPECT_EQ(
            subdomain,
            ProductMDomain(mesh_x_vx, MCoord<MeshX, MeshVx>(1, 0), MLength<MeshX, MeshVx>(1, 4)));
}
