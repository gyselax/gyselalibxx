#include <gtest/gtest.h>

#include "bsplines_non_uniform.h"

using NonUniformMeshX = NonUniformMesh<Dim::X>;

TEST(BSplinesNonUniform, Constructor)
{
    std::vector<double> breaks({0.0, 0.5, 1.0, 1.5, 2.0});

    NonUniformMeshX mesh(breaks, MCoordX(0));
    MDomainImpl<NonUniformMeshX> const dom(mesh, MCoordX(5));
    NonUniformBSplines bsplines(2, dom);
    EXPECT_EQ(bsplines.degree(), 2);
    EXPECT_EQ(bsplines.is_periodic(), Dim::X::PERIODIC);
    EXPECT_EQ(bsplines.xmin(), 0.);
    EXPECT_EQ(bsplines.xmax(), 2.);
    EXPECT_EQ(bsplines.ncells(), 4);
}
