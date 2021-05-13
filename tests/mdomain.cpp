#include <gtest/gtest.h>

#include "mdomain.h"

TEST(MDomain, RangeFor)
{
    MDomainX dom(0., 10., 0, 2);
    int ii = 0;
    for (auto&& x : dom) {
        ASSERT_LE(0, x);
        ASSERT_EQ(x, ii);
        ASSERT_LT(x, 2);
        ++ii;
    }
}

TEST(MDomain, lbound)
{
    MDomainXVx const
            dom2d(RCoordXVx(0., 0.), RCoordXVx(1., 1.), MCoordXVx(0, 0), MCoordXVx(100, 100));
    ASSERT_EQ(dom2d.to_real(dom2d.lbound()).get<Dim::X>(), 0.);
}
