#include <iosfwd>
#include <memory>

#include <ddc/mcoord.h>
#include <ddc/non_uniform_mesh.h>
#include <ddc/product_mdomain.h>
#include <ddc/product_mesh.h>
#include <ddc/rcoord.h>
#include <ddc/taggedvector.h>
#include <ddc/uniform_mesh.h>

#include <gtest/gtest.h>

class DimX;
class DimVx;

using MeshX = UniformMesh<DimX>;
using MeshVx = NonUniformMesh<DimVx>;

class ProductMDomainTest : public ::testing::Test
{
protected:
    std::size_t npoints = 11;
    MeshX mesh_x = MeshX(2., 3., npoints);
    MeshVx mesh_vx = MeshVx({-1., 0., 2., 4.});
    ProductMesh<MeshX, MeshVx> mesh_x_vx = ProductMesh(mesh_x, mesh_vx);
    ProductMDomain<MeshX, MeshVx> domain_x_vx = ProductMDomain(
            mesh_x_vx,
            MCoord<MeshX, MeshVx>(0, 0),
            MCoord<MeshX, MeshVx>(npoints - 1, 3));
};

TEST_F(ProductMDomainTest, constructor)
{
    EXPECT_EQ(domain_x_vx.extents(), (MCoord<MeshX, MeshVx>(11, 4)));
    EXPECT_EQ(domain_x_vx.front(), (MCoord<MeshX, MeshVx>(0, 0)));
    EXPECT_EQ(domain_x_vx.back(), (MCoord<MeshX, MeshVx>(10, 3)));
    EXPECT_EQ(domain_x_vx.size(), 11 * 4);
}

TEST_F(ProductMDomainTest, rmin_rmax)
{
    EXPECT_EQ(domain_x_vx.rmin(), (RCoord<DimX, DimVx>(2., -1.)));
    EXPECT_EQ(domain_x_vx.rmax(), (RCoord<DimX, DimVx>(3., 4.)));
}

TEST_F(ProductMDomainTest, subdomain)
{
    ProductMDomain subdomain_x(ProductMesh<MeshX>(mesh_x), MCoord<MeshX>(1), MCoord<MeshX>(1));
    auto subdomain = domain_x_vx.restrict(subdomain_x);
    EXPECT_EQ(
            subdomain,
            ProductMDomain(mesh_x_vx, MCoord<MeshX, MeshVx>(1, 0), MCoord<MeshX, MeshVx>(1, 3)));
}
