#include <gtest/gtest.h>

#include "mdomain.h"

TEST(MDomain, RangeFor)
{
    MDomainX dom(0., 10., 0ul, 2ul);
    int ii = 0;
    for (auto&& x : dom) {
        ASSERT_LE(0, x);
        EXPECT_EQ(x, ii);
        ASSERT_LT(x, 2);
        ++ii;
    }
}

TEST(MDomain, ubound)
{
    MDomainXVx const
            dom2d(RCoordXVx(0., 0.),
                  RCoordXVx(2., 2.),
                  MCoordXVx(0ul, 0ul),
                  MCoordXVx(100ul, 200ul));
    EXPECT_EQ(dom2d.ubound().get<Dim::X>(), 100ul);
    EXPECT_EQ(dom2d.ubound<Dim::X>(), 100ul);

    EXPECT_EQ(dom2d.ubound().get<Dim::Vx>(), 200ul);
    EXPECT_EQ(dom2d.ubound<Dim::Vx>(), 200ul);
}

TEST(MDomain, rmax)
{
    MDomainXVx const
            dom2d(RCoordXVx(0., 0.),
                  RCoordXVx(2., 2.),
                  MCoordXVx(0ul, 0ul),
                  MCoordXVx(100ul, 200ul));
    EXPECT_EQ(dom2d.mesh().to_real(dom2d.ubound()).get<Dim::X>(), 200.);
    EXPECT_EQ(dom2d.rmax().get<Dim::X>(), 200.);
    EXPECT_EQ(dom2d.rmax<Dim::X>(), 200.);

    EXPECT_EQ(dom2d.mesh().to_real(dom2d.ubound()).get<Dim::Vx>(), 400.);
    EXPECT_EQ(dom2d.rmax().get<Dim::Vx>(), 400.);
    EXPECT_EQ(dom2d.rmax<Dim::Vx>(), 400.);
}
