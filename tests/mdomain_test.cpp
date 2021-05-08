#include <gtest/gtest.h>

#include "mdomain.h"

TEST(MDomain, RangeFor)
{
    MDomainX dom(0., 10., 0, 2);
    int ii=0;
    for ( auto&& x: dom) {
        ASSERT_LE(0, x);
        ASSERT_EQ(x, ii);
        ASSERT_LT(x, 2);
        ++ii;
    }
}
