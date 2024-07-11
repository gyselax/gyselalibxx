// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace DDC_HIP_5_7_ANONYMOUS_NAMESPACE_WORKAROUND(DISCRETE_DOMAIN_CPP)
{
    struct DDimX
    {
    };
    using DElemX = ddc::DiscreteElement<DDimX>;
    using DVectX = ddc::DiscreteVector<DDimX>;
    using DDomX = ddc::DiscreteDomain<DDimX>;


    struct DDimY
    {
    };
    using DElemY = ddc::DiscreteElement<DDimY>;
    using DVectY = ddc::DiscreteVector<DDimY>;
    using DDomY = ddc::DiscreteDomain<DDimY>;


    struct DDimZ
    {
    };
    using DElemZ = ddc::DiscreteElement<DDimZ>;
    using DVectZ = ddc::DiscreteVector<DDimZ>;
    using DDomZ = ddc::DiscreteDomain<DDimZ>;


    using DElemXY = ddc::DiscreteElement<DDimX, DDimY>;
    using DVectXY = ddc::DiscreteVector<DDimX, DDimY>;
    using DDomXY = ddc::DiscreteDomain<DDimX, DDimY>;


    using DElemYX = ddc::DiscreteElement<DDimY, DDimX>;
    using DVectYX = ddc::DiscreteVector<DDimY, DDimX>;
    using DDomYX = ddc::DiscreteDomain<DDimY, DDimX>;

    using DElemXZ = ddc::DiscreteElement<DDimX, DDimZ>;
    using DVectXZ = ddc::DiscreteVector<DDimX, DDimZ>;
    using DDomXZ = ddc::DiscreteDomain<DDimX, DDimZ>;

    using DElemZY = ddc::DiscreteElement<DDimZ, DDimY>;
    using DVectZY = ddc::DiscreteVector<DDimZ, DDimY>;
    using DDomZY = ddc::DiscreteDomain<DDimZ, DDimY>;


    using DElemXYZ = ddc::DiscreteElement<DDimX, DDimY, DDimZ>;
    using DVectXYZ = ddc::DiscreteVector<DDimX, DDimY, DDimZ>;
    using DDomXYZ = ddc::DiscreteDomain<DDimX, DDimY, DDimZ>;

    using DElemZYX = ddc::DiscreteElement<DDimZ, DDimY, DDimX>;
    using DVectZYX = ddc::DiscreteVector<DDimZ, DDimY, DDimX>;
    using DDomZYX = ddc::DiscreteDomain<DDimZ, DDimY, DDimX>;

    static DElemX constexpr lbound_x(50);
    static DVectX constexpr nelems_x(3);
    static DElemX constexpr sentinel_x(lbound_x + nelems_x);
    static DElemX constexpr ubound_x(sentinel_x - 1);


    static DElemY constexpr lbound_y(4);
    static DVectY constexpr nelems_y(12);
    static DElemY constexpr sentinel_y(lbound_y + nelems_y);
    static DElemY constexpr ubound_y(sentinel_y - 1);

    static DElemZ constexpr lbound_z(7);
    static DVectZ constexpr nelems_z(15);

    static DElemXY constexpr lbound_x_y(lbound_x, lbound_y);
    static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);
    static DElemXY constexpr ubound_x_y(ubound_x, ubound_y);

    static DElemXZ constexpr lbound_x_z(lbound_x, lbound_z);
    static DVectXZ constexpr nelems_x_z(nelems_x, nelems_z);

} // namespace )

TEST(DiscreteDomainTest, Constructor)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(dom_x_y.extents(), nelems_x_y);
    EXPECT_EQ(dom_x_y.front(), lbound_x_y);
    EXPECT_EQ(dom_x_y.back(), ubound_x_y);
    EXPECT_EQ(dom_x_y.size(), nelems_x.value() * nelems_y.value());

    DDomX const dom_x(lbound_x, nelems_x);
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

TEST(DiscreteDomainTest, ConstructorFromDiscreteDomains)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    DDomZ const dom_z(lbound_z, nelems_z);
    DDomXYZ const dom_x_y_z(dom_z, dom_x_y);
    EXPECT_EQ(dom_x_y_z.front(), DElemXYZ(lbound_x, lbound_y, lbound_z));
    EXPECT_EQ(dom_x_y_z.extents(), DVectXYZ(nelems_x, nelems_y, nelems_z));
}

TEST(DiscreteDomainTest, EmptyDomain)
{
    DDomXY const dom_x_y = DDomXY();
    EXPECT_EQ(dom_x_y.extents(), DVectXY(0, 0));
    EXPECT_EQ(dom_x_y.size(), 0);
    EXPECT_TRUE(dom_x_y.empty());
}

TEST(DiscreteDomainTest, CompareSameDomains)
{
    DDomXY const dom_x_y_1(lbound_x_y, nelems_x_y);
    DDomXY const dom_x_y_2(dom_x_y_1);
    EXPECT_TRUE(dom_x_y_1 == dom_x_y_2);
    EXPECT_TRUE(dom_x_y_1 == DDomYX(dom_x_y_2));
    EXPECT_FALSE(dom_x_y_1 != dom_x_y_2);
    EXPECT_FALSE(dom_x_y_1 != DDomYX(dom_x_y_2));
}

TEST(DiscreteDomainTest, CompareDifferentDomains)
{
    DDomXY const dom_x_y_1(DElemXY(0, 1), DVectXY(1, 2));
    DDomXY const dom_x_y_2(DElemXY(2, 3), DVectXY(3, 4));
    EXPECT_FALSE(dom_x_y_1 == dom_x_y_2);
    EXPECT_FALSE(dom_x_y_1 == DDomYX(dom_x_y_2));
    EXPECT_TRUE(dom_x_y_1 != dom_x_y_2);
    EXPECT_TRUE(dom_x_y_1 != DDomYX(dom_x_y_2));
}

TEST(DiscreteDomainTest, CompareEmptyDomains)
{
    DDomXY const dom_x_y_1(DElemXY(4, 1), DVectXY(0, 0));
    DDomXY const dom_x_y_2(DElemXY(3, 9), DVectXY(0, 0));
    EXPECT_TRUE(dom_x_y_1.empty());
    EXPECT_TRUE(dom_x_y_2.empty());
    EXPECT_TRUE(dom_x_y_1 == dom_x_y_2);
    EXPECT_FALSE(dom_x_y_1 != dom_x_y_2);
}

TEST(DiscreteDomainTest, Subdomain)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    ddc::DiscreteElement<DDimX> const lbound_subdomain_x(lbound_x + 1);
    ddc::DiscreteVector<DDimX> const npoints_subdomain_x(nelems_x - 2);
    DDomX const subdomain_x(lbound_subdomain_x, npoints_subdomain_x);
    DDomXY const subdomain = dom_x_y.restrict(subdomain_x);
    EXPECT_EQ(
            subdomain,
            DDomXY(ddc::DiscreteElement<DDimX, DDimY>(lbound_subdomain_x, lbound_y),
                   ddc::DiscreteVector<DDimX, DDimY>(npoints_subdomain_x, nelems_y)));
}

TEST(DiscreteDomainTest, RangeFor)
{
    DDomX const dom(lbound_x, nelems_x);
    DElemX ii = lbound_x;
    for (DElemX ix : dom) {
        EXPECT_LE(lbound_x, ix);
        EXPECT_EQ(ix, ii);
        EXPECT_LE(ix, ubound_x);
        ++ii.uid<DDimX>();
    }
}

TEST(DiscreteDomainTest, DiffEmpty)
{
    DDomX const dom_x = DDomX();
    auto const subdomain1 = ddc::remove_dims_of(dom_x, dom_x);
    auto const subdomain2 = ddc::remove_dims_of<DDimX>(dom_x);
    EXPECT_EQ(subdomain1, ddc::DiscreteDomain<>());
    EXPECT_EQ(subdomain2, ddc::DiscreteDomain<>());
}

TEST(DiscreteDomainTest, Diff)
{
    DDomX const dom_x = DDomX();
    DDomXY const dom_x_y = DDomXY();
    DDomZY const dom_z_y = DDomZY();
    auto const subdomain1 = ddc::remove_dims_of(dom_x_y, dom_z_y);
    auto const subdomain2 = ddc::remove_dims_of<DDimZ, DDimY>(dom_x_y);
    EXPECT_EQ(subdomain1, dom_x);
    EXPECT_EQ(subdomain2, dom_x);
}

TEST(DiscreteDomainTest, Replace)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    DDomZ const dom_z(lbound_z, nelems_z);
    DDomXZ const dom_x_z(lbound_x_z, nelems_x_z);
    auto const subdomain = ddc::replace_dim_of<DDimY, DDimZ>(dom_x_y, dom_z);
    EXPECT_EQ(subdomain, dom_x_z);
}


TEST(DiscreteDomainTest, TakeFirst)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(dom_x_y.take_first(DVectXY(1, 4)), DDomXY(dom_x_y.front(), DVectXY(1, 4)));
}

TEST(DiscreteDomainTest, TakeLast)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.take_last(DVectXY(1, 4)),
            DDomXY(dom_x_y.front() + dom_x_y.extents() - DVectXY(1, 4), DVectXY(1, 4)));
}

TEST(DiscreteDomainTest, RemoveFirst)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.remove_first(DVectXY(1, 4)),
            DDomXY(dom_x_y.front() + DVectXY(1, 4), dom_x_y.extents() - DVectXY(1, 4)));
}

TEST(DiscreteDomainTest, RemoveLast)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.remove_last(DVectXY(1, 4)),
            DDomXY(dom_x_y.front(), dom_x_y.extents() - DVectXY(1, 4)));
}

TEST(DiscreteDomainTest, Remove)
{
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    EXPECT_EQ(
            dom_x_y.remove(DVectXY(1, 4), DVectXY(1, 1)),
            DDomXY(dom_x_y.front() + DVectXY(1, 4), dom_x_y.extents() - DVectXY(2, 5)));
}

TEST(DiscreteDomainTest, SliceDomainXTooearly)
{
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    DDomX const subdomain_x(lbound_x - 1, nelems_x);
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            dom_x_y.restrict(subdomain_x),
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_begin\).*uid<ODDims>\(odomain\.m_element_begin\))rgx");
#else
    GTEST_SKIP();
#endif
}

TEST(DiscreteDomainTest, SliceDomainXToolate)
{
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    DDomX const subdomain_x(lbound_x, nelems_x + 1);
    DDomXY const dom_x_y(lbound_x_y, nelems_x_y);
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            dom_x_y.restrict(subdomain_x),
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_end\).*uid<ODDims>\(odomain\.m_element_end\).*)rgx");
#else
    GTEST_SKIP();
#endif
}

TEST(DiscreteDomainTest, Transpose3DConstructor)
{
    DDomX const dom_x(lbound_x, nelems_x);
    DDomY const dom_y(lbound_y, nelems_y);
    DDomZ const dom_z(lbound_z, nelems_z);
    DDomXYZ const dom_x_y_z(dom_x, dom_y, dom_z);
    DDomZYX const dom_z_y_x(dom_x_y_z);
    EXPECT_EQ(ddc::select<DDimX>(dom_x_y_z.front()), ddc::select<DDimX>(dom_z_y_x.front()));
    EXPECT_EQ(ddc::select<DDimY>(dom_x_y_z.front()), ddc::select<DDimY>(dom_z_y_x.front()));
    EXPECT_EQ(ddc::select<DDimZ>(dom_x_y_z.front()), ddc::select<DDimZ>(dom_z_y_x.front()));
    EXPECT_EQ(ddc::select<DDimX>(dom_x_y_z.back()), ddc::select<DDimX>(dom_z_y_x.back()));
    EXPECT_EQ(ddc::select<DDimY>(dom_x_y_z.back()), ddc::select<DDimY>(dom_z_y_x.back()));
    EXPECT_EQ(ddc::select<DDimZ>(dom_x_y_z.back()), ddc::select<DDimZ>(dom_z_y_x.back()));
}
