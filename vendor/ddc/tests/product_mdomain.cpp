// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>

#include <gtest/gtest.h>

class DimX;
using RCoordX = RCoord<DimX>;

class DimVx;
using RCoordVx = RCoord<DimVx>;

using MeshX = UniformMesh<DimX>;
using MCoordX = MCoord<MeshX>;
using MLengthX = MLength<MeshX>;

using MDomainX = ProductMDomain<MeshX>;

using MeshVx = NonUniformMesh<DimVx>;
using MCoordVx = MCoord<MeshVx>;
using MLengthVx = MLength<MeshVx>;

using MCoordXVx = MCoord<MeshX, MeshVx>;
using MLengthXVx = MLength<MeshX, MeshVx>;

using MDomainSpXVx = ProductMDomain<MeshX, MeshVx>;


RCoordX constexpr origin_x(0);
RCoordX constexpr step_x(.01);
MeshX constexpr mesh_x = MeshX(origin_x, step_x);

MCoordX constexpr lbound_x(50);
MLengthX constexpr npoints_x(51);
MCoordX constexpr ubound_x(lbound_x + npoints_x - 1);

MLengthVx constexpr npoints_vx(4);
std::array<RCoordVx, npoints_vx> constexpr coords_vx
        = {RCoordVx(-1.), RCoordVx(0.), RCoordVx(2.), RCoordVx(4.)};
MCoordVx constexpr lbound_vx(0);
MeshVx const mesh_vx = MeshVx(coords_vx);

MCoordXVx constexpr lbound_x_vx = MCoord<MeshX, MeshVx> {lbound_x, lbound_vx};
MLengthXVx constexpr npoints_x_vx = MLength<MeshX, MeshVx>(npoints_x, npoints_vx);
MLengthXVx constexpr ubound_x_vx
        = MLength<MeshX, MeshVx>(lbound_x + npoints_x - 1, lbound_vx + npoints_vx - 1);


TEST(ProductMDomainTest, constructor)
{
    MDomainSpXVx const dom_x_vx = MDomainSpXVx(mesh_x, mesh_vx, lbound_x_vx, npoints_x_vx);
    EXPECT_EQ(dom_x_vx.extents(), npoints_x_vx);
    EXPECT_EQ(dom_x_vx.front(), lbound_x_vx);
    EXPECT_EQ(dom_x_vx.back(), ubound_x_vx);
    EXPECT_EQ(dom_x_vx.size(), npoints_x.value() * npoints_vx.value());

    MDomainX const dom_x = MDomainX(mesh_x, lbound_x, npoints_x);
    EXPECT_EQ(dom_x.mesh<MeshX>(), mesh_x);
    EXPECT_EQ(dom_x.size(), npoints_x);
    EXPECT_EQ(dom_x.empty(), false);
    EXPECT_EQ(dom_x[0], lbound_x);
    EXPECT_EQ(dom_x.front(), lbound_x);
    EXPECT_EQ(dom_x.back(), ubound_x);

    MDomainX const empty_domain(mesh_x, lbound_x, MLengthX(0));
    EXPECT_EQ(empty_domain.mesh<MeshX>(), mesh_x);
    EXPECT_EQ(empty_domain.size(), 0);
    EXPECT_EQ(empty_domain.empty(), true);
    EXPECT_EQ(empty_domain[0], lbound_x);
}

TEST(ProductMDomainTest, rmin_rmax)
{
    MDomainX const dom_x(mesh_x, lbound_x, npoints_x);
    EXPECT_EQ(dom_x.rmin(), lbound_x.value() * step_x + origin_x);
    EXPECT_EQ(dom_x.rmax(), ubound_x.value() * step_x + origin_x);

    MDomainSpXVx const dom_x_vx(mesh_x, mesh_vx, lbound_x_vx, npoints_x_vx);
    RCoord<DimX, DimVx> const
            rmin_vx(RCoordX(lbound_x.value() * step_x) + origin_x, coords_vx[lbound_vx]);
    RCoord<DimX, DimVx> const
            rmax_vx(RCoordX(ubound_x.value() * step_x) + origin_x, coords_vx[npoints_vx - 1]);
    EXPECT_EQ(dom_x_vx.rmin(), rmin_vx);
    EXPECT_EQ(dom_x_vx.rmax(), rmax_vx);
}

TEST(ProductMDomainTest, subdomain)
{
    MDomainSpXVx const dom_x_vx = MDomainSpXVx(mesh_x, mesh_vx, lbound_x_vx, npoints_x_vx);
    MCoord<MeshX> const lbound_subdomain_x(lbound_x + 1);
    MLength<MeshX> const npoints_subdomain_x(npoints_x - 2);
    MDomainX const subdomain_x(mesh_x, lbound_subdomain_x, npoints_subdomain_x);
    MDomainSpXVx const subdomain = dom_x_vx.restrict(subdomain_x);
    EXPECT_EQ(
            subdomain,
            MDomainSpXVx(
                    mesh_x,
                    mesh_vx,
                    MCoord<MeshX, MeshVx>(lbound_subdomain_x, lbound_vx),
                    MLength<MeshX, MeshVx>(npoints_subdomain_x, npoints_vx)));
}

TEST(ProductMDomainTest, RangeFor)
{
    MDomainX const dom = MDomainX(mesh_x, lbound_x, npoints_x);
    MCoordX ii = lbound_x;
    for (MCoordX ix : dom) {
        ASSERT_LE(lbound_x, ix);
        EXPECT_EQ(ix, ii);
        ASSERT_LE(ix, ubound_x);
        ++ii.get<MeshX>();
    }
}
