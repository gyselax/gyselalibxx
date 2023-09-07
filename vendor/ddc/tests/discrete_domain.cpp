// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

struct DDimX;
using DElemX = ddc::DiscreteElement<DDimX>;
using DVectX = ddc::DiscreteVector<DDimX>;
using DDomX = ddc::DiscreteDomain<DDimX>;


struct DDimY;
using DElemY = ddc::DiscreteElement<DDimY>;
using DVectY = ddc::DiscreteVector<DDimY>;
using DDomY = ddc::DiscreteDomain<DDimY>;


struct DDimZ;
using DElemZ = ddc::DiscreteElement<DDimZ>;
using DVectZ = ddc::DiscreteVector<DDimZ>;
using DDomZ = ddc::DiscreteDomain<DDimZ>;


using DElemXY = ddc::DiscreteElement<DDimX, DDimY>;
using DVectXY = ddc::DiscreteVector<DDimX, DDimY>;
using DDomXY = ddc::DiscreteDomain<DDimX, DDimY>;


using DElemYX = ddc::DiscreteElement<DDimY, DDimX>;
using DVectYX = ddc::DiscreteVector<DDimY, DDimX>;
using DDomYX = ddc::DiscreteDomain<DDimY, DDimX>;

using DElemZY = ddc::DiscreteElement<DDimZ, DDimY>;
using DVectZY = ddc::DiscreteVector<DDimZ, DDimY>;
using DDomZY = ddc::DiscreteDomain<DDimZ, DDimY>;



static DElemX constexpr lbound_x(50);
static DVectX constexpr nelems_x(3);
static DElemX constexpr sentinel_x(lbound_x + nelems_x);
static DElemX constexpr ubound_x(sentinel_x - 1); //TODO: correct type


static DElemY constexpr lbound_y(4);
static DVectY constexpr nelems_y(12);
static DElemY constexpr sentinel_y(lbound_y + nelems_y);
static DElemY constexpr ubound_y(sentinel_y - 1); //TODO: correct type


static DElemXY constexpr lbound_x_y(lbound_x, lbound_y);
static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);
static DElemXY constexpr ubound_x_y(ubound_x, ubound_y);

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

TEST(ProductMDomainTest, EmptyDomain)
{
    DDomXY const dom_x_y = DDomXY();
    EXPECT_EQ(dom_x_y.extents(), DVectXY(0, 0));
    EXPECT_EQ(dom_x_y.size(), 0);
    EXPECT_TRUE(dom_x_y.empty());
}

TEST(ProductMDomainTest, Subdomain)
{
    DDomXY const dom_x_y = DDomXY(lbound_x_y, nelems_x_y);
    ddc::DiscreteElement<DDimX> const lbound_subdomain_x(lbound_x + 1);
    ddc::DiscreteVector<DDimX> const npoints_subdomain_x(nelems_x - 2);
    DDomX const subdomain_x(lbound_subdomain_x, npoints_subdomain_x);
    DDomXY const subdomain = dom_x_y.restrict(subdomain_x);
    EXPECT_EQ(
            subdomain,
            DDomXY(ddc::DiscreteElement<DDimX, DDimY>(lbound_subdomain_x, lbound_y),
                   ddc::DiscreteVector<DDimX, DDimY>(npoints_subdomain_x, nelems_y)));
}

TEST(ProductMDomainTest, RangeFor)
{
    DDomX const dom = DDomX(lbound_x, nelems_x);
    DElemX ii = lbound_x;
    for (DElemX ix : dom) {
        EXPECT_LE(lbound_x, ix);
        EXPECT_EQ(ix, ii);
        EXPECT_LE(ix, ubound_x);
        ++ii.uid<DDimX>();
    }
}

TEST(ProductMDomainTest, DiffEmpty)
{
    DDomX const dom_x = DDomX();
    auto const subdomain = ddc::remove_dims_of(dom_x, dom_x);
    EXPECT_EQ(subdomain, ddc::DiscreteDomain<>());
}

TEST(ProductMDomainTest, Diff)
{
    DDomX const dom_x = DDomX();
    DDomXY const dom_x_y = DDomXY();
    DDomZY const dom_z_y = DDomZY();
    auto const subdomain = ddc::remove_dims_of(dom_x_y, dom_z_y);
    EXPECT_EQ(subdomain, dom_x);
}

TEST(ProductMDomainTest, TakeFirst)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(dom_x_y.take_first(DVectXY(1, 4)), DDomXY(dom_x_y.front(), DVectXY(1, 4)));
}

TEST(ProductMDomainTest, TakeLast)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.take_last(DVectXY(1, 4)),
            DDomXY(dom_x_y.front() + dom_x_y.extents() - DVectXY(1, 4), DVectXY(1, 4)));
}

TEST(ProductMDomainTest, RemoveFirst)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.remove_first(DVectXY(1, 4)),
            DDomXY(dom_x_y.front() + DVectXY(1, 4), dom_x_y.extents() - DVectXY(1, 4)));
}

TEST(ProductMDomainTest, RemoveLast)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.remove_last(DVectXY(1, 4)),
            DDomXY(dom_x_y.front(), dom_x_y.extents() - DVectXY(1, 4)));
}

TEST(ProductMDomainTest, Remove)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.remove(DVectXY(1, 4), DVectXY(1, 1)),
            DDomXY(dom_x_y.front() + DVectXY(1, 4), dom_x_y.extents() - DVectXY(2, 5)));
}
