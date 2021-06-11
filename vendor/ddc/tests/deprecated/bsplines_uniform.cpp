#include <gtest/gtest.h>

#include "deprecated/bsplines_uniform.h"

struct DimX
{
    static constexpr bool PERIODIC = true;
};
using UniformMDomainX = UniformMDomain<DimX>;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<DimX>;

namespace deprecated {

TEST(BSplinesUniform, Constructor)
{
    UniformMDomainX const dom(RCoordX(0.), RCoordX(2.02), MCoordX(0), MCoordX(101));
    UniformBSplines bsplines(2, dom);
    EXPECT_EQ(bsplines.degree(), 2);
    EXPECT_EQ(bsplines.is_periodic(), DimX::PERIODIC);
    EXPECT_EQ(bsplines.xmin(), 0.);
    EXPECT_EQ(bsplines.xmax(), 2.);
    EXPECT_EQ(bsplines.ncells(), 100);
}

} // namespace deprecated
