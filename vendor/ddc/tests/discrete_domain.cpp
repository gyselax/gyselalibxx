// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

struct DDimX;
using ElemX = DiscreteCoordinate<DDimX>;
using DVectX = DiscreteVector<DDimX>;
using DDomX = DiscreteDomain<DDimX>;


struct DDimY;
using ElemY = DiscreteCoordinate<DDimY>;
using DVectY = DiscreteVector<DDimY>;
using DDomY = DiscreteDomain<DDimY>;


struct DDimZ;
using ElemZ = DiscreteCoordinate<DDimZ>;
using DVectZ = DiscreteVector<DDimZ>;
using DDomZ = DiscreteDomain<DDimZ>;


using ElemXY = DiscreteCoordinate<DDimX, DDimY>;
using DVectXY = DiscreteVector<DDimX, DDimY>;
using DDomXY = DiscreteDomain<DDimX, DDimY>;


using ElemYX = DiscreteCoordinate<DDimY, DDimX>;
using DVectYX = DiscreteVector<DDimY, DDimX>;
using DDomYX = DiscreteDomain<DDimY, DDimX>;


static ElemX constexpr lbound_x(50);
static DVectX constexpr nelems_x(3);
static ElemX constexpr sentinel_x = lbound_x + nelems_x;
static ElemX constexpr ubound_x = ElemX(sentinel_x - 1); //TODO: correct type


static ElemY constexpr lbound_y(4);
static DVectY constexpr nelems_y(12);
static ElemY constexpr sentinel_y = lbound_y + nelems_y;
static ElemY constexpr ubound_y = ElemY(sentinel_y - 1); //TODO: correct type


static ElemXY constexpr lbound_x_y {lbound_x, lbound_y};
static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);
static ElemXY constexpr ubound_x_y(ubound_x, ubound_y);

} // namespace

TEST(ProductMDomainTest, Constructor)
{
    DDomXY const dom_x_y = DDomXY(lbound_x_y, nelems_x_y);
    EXPECT_EQ(dom_x_y.extents(), nelems_x_y);
    EXPECT_EQ(dom_x_y.front(), lbound_x_y);
    EXPECT_EQ(dom_x_y.back(), ubound_x_y);
    EXPECT_EQ(dom_x_y.size(), nelems_x.value() * nelems_y.value());

    DDomX const dom_x = DDomX(lbound_x, nelems_x);
    EXPECT_EQ(dom_x.size(), nelems_x);
    EXPECT_EQ(dom_x.empty(), false);
    EXPECT_EQ(dom_x[0], lbound_x);
    EXPECT_EQ(dom_x.front(), lbound_x);
    EXPECT_EQ(dom_x.back(), ubound_x);

    DDomX const empty_domain(lbound_x, DVectX(0));
    EXPECT_EQ(empty_domain.size(), 0);
    EXPECT_EQ(empty_domain.empty(), true);
    EXPECT_EQ(empty_domain[0], lbound_x);
}

TEST(ProductMDomainTest, Subdomain)
{
    DDomXY const dom_x_y = DDomXY(lbound_x_y, nelems_x_y);
    DiscreteCoordinate<DDimX> const lbound_subdomain_x(lbound_x + 1);
    DiscreteVector<DDimX> const npoints_subdomain_x(nelems_x - 2);
    DDomX const subdomain_x(lbound_subdomain_x, npoints_subdomain_x);
    DDomXY const subdomain = dom_x_y.restrict(subdomain_x);
    EXPECT_EQ(
            subdomain,
            DDomXY(DiscreteCoordinate<DDimX, DDimY>(lbound_subdomain_x, lbound_y),
                   DiscreteVector<DDimX, DDimY>(npoints_subdomain_x, nelems_y)));
}

TEST(ProductMDomainTest, RangeFor)
{
    DDomX const dom = DDomX(lbound_x, nelems_x);
    ElemX ii = lbound_x;
    for (ElemX ix : dom) {
        ASSERT_LE(lbound_x, ix);
        EXPECT_EQ(ix, ii);
        ASSERT_LE(ix, ubound_x);
        ++ii.get<DDimX>();
    }
}
