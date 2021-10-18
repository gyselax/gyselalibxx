// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/NonUniformDiscretization>
#include <ddc/UniformDiscretization>

#include <gtest/gtest.h>

namespace {

class RDimX;
using CoordX = Coordinate<RDimX>;

class DimY;
using RCoordY = Coordinate<DimY>;

using IDimX = UniformDiscretization<RDimX>;
using IndexX = DiscreteCoordinate<IDimX>;
using IVectX = DiscreteVector<IDimX>;

using IDomainX = DiscreteDomain<IDimX>;

using DDimY = NonUniformDiscretization<DimY>;
using MCoordY = DiscreteCoordinate<DDimY>;
using MLengthY = DiscreteVector<DDimY>;

using MCoordXY = DiscreteCoordinate<IDimX, DDimY>;
using MLengthXY = DiscreteVector<IDimX, DDimY>;

using MDomainXY = DiscreteDomain<IDimX, DDimY>;


CoordX constexpr origin_x(0);
CoordX constexpr step_x(.01);
IDimX constexpr idim_x = IDimX(origin_x, step_x);

IndexX constexpr lbound_x(50);
IVectX constexpr npoints_x(51);
IndexX constexpr ubound_x(lbound_x + npoints_x - 1);

MLengthY constexpr npoints_y(4);
std::array<RCoordY, npoints_y> constexpr coords_y
        = {RCoordY(-1.), RCoordY(0.), RCoordY(2.), RCoordY(4.)};
MCoordY constexpr lbound_y(0);
DDimY const ddim_y = DDimY(coords_y);

MCoordXY constexpr lbound_x_y {lbound_x, lbound_y};
MLengthXY constexpr npoints_x_y(npoints_x, npoints_y);
MLengthXY constexpr ubound_x_y(lbound_x + npoints_x - 1, lbound_y + npoints_y - 1);

} // namespace

TEST(ProductMDomainTest, constructor)
{
    MDomainXY const dom_x_y = MDomainXY(idim_x, ddim_y, lbound_x_y, npoints_x_y);
    EXPECT_EQ(dom_x_y.extents(), npoints_x_y);
    EXPECT_EQ(dom_x_y.front(), lbound_x_y);
    EXPECT_EQ(dom_x_y.back(), ubound_x_y);
    EXPECT_EQ(dom_x_y.size(), npoints_x.value() * npoints_y.value());

    IDomainX const dom_x = IDomainX(idim_x, lbound_x, npoints_x);
    EXPECT_EQ(dom_x.mesh<IDimX>(), idim_x);
    EXPECT_EQ(dom_x.size(), npoints_x);
    EXPECT_EQ(dom_x.empty(), false);
    EXPECT_EQ(dom_x[0], lbound_x);
    EXPECT_EQ(dom_x.front(), lbound_x);
    EXPECT_EQ(dom_x.back(), ubound_x);

    IDomainX const empty_domain(idim_x, lbound_x, IVectX(0));
    EXPECT_EQ(empty_domain.mesh<IDimX>(), idim_x);
    EXPECT_EQ(empty_domain.size(), 0);
    EXPECT_EQ(empty_domain.empty(), true);
    EXPECT_EQ(empty_domain[0], lbound_x);
}

TEST(ProductMDomainTest, rmin_rmax)
{
    IDomainX const dom_x(idim_x, lbound_x, npoints_x);
    EXPECT_EQ(dom_x.rmin(), lbound_x.value() * step_x + origin_x);
    EXPECT_EQ(dom_x.rmax(), ubound_x.value() * step_x + origin_x);

    MDomainXY const dom_x_y(idim_x, ddim_y, lbound_x_y, npoints_x_y);
    Coordinate<RDimX, DimY> const
            rmin_y(CoordX(lbound_x.value() * step_x) + origin_x, coords_y[lbound_y]);
    Coordinate<RDimX, DimY> const
            rmax_y(CoordX(ubound_x.value() * step_x) + origin_x, coords_y[npoints_y - 1]);
    EXPECT_EQ(dom_x_y.rmin(), rmin_y);
    EXPECT_EQ(dom_x_y.rmax(), rmax_y);
}

TEST(ProductMDomainTest, subdomain)
{
    MDomainXY const dom_x_y = MDomainXY(idim_x, ddim_y, lbound_x_y, npoints_x_y);
    DiscreteCoordinate<IDimX> const lbound_subdomain_x(lbound_x + 1);
    DiscreteVector<IDimX> const npoints_subdomain_x(npoints_x - 2);
    IDomainX const subdomain_x(idim_x, lbound_subdomain_x, npoints_subdomain_x);
    MDomainXY const subdomain = dom_x_y.restrict(subdomain_x);
    EXPECT_EQ(
            subdomain,
            MDomainXY(
                    idim_x,
                    ddim_y,
                    DiscreteCoordinate<IDimX, DDimY>(lbound_subdomain_x, lbound_y),
                    DiscreteVector<IDimX, DDimY>(npoints_subdomain_x, npoints_y)));
}

TEST(ProductMDomainTest, RangeFor)
{
    IDomainX const dom = IDomainX(idim_x, lbound_x, npoints_x);
    IndexX ii = lbound_x;
    for (IndexX ix : dom) {
        ASSERT_LE(lbound_x, ix);
        EXPECT_EQ(ix, ii);
        ASSERT_LE(ix, ubound_x);
        ++ii.get<IDimX>();
    }
}
