#include <gtest/gtest.h>

#include "bsplines_uniform.h"

TEST(BSplinesUniform, Constructor)
{
    UniformMDomainX const dom(RCoordX(0.), RCoordX(2.02), MCoordX(0), MCoordX(101));
    UniformBSplines bsplines(2, dom);
    EXPECT_EQ(bsplines.degree(), 2);
    EXPECT_EQ(bsplines.is_periodic(), Dim::X::PERIODIC);
    EXPECT_EQ(bsplines.xmin(), 0.);
    EXPECT_EQ(bsplines.xmax(), 2.);
    EXPECT_EQ(bsplines.ncells(), 100);
}
