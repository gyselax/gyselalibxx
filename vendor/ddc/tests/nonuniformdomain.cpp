#include <gtest/gtest.h>

#include "mdomain.h"
#include "non_uniform_mesh.h"

class DimX;
using NonUniformMeshX = NonUniformMesh<DimX>;
using NonUniformDomainX = MDomainImpl<NonUniformMeshX>;
using MCoordX = MCoord<DimX>;

TEST(NonUniformDomainXTest, Constructor)
{
    std::vector<double> breaks(4);
    for (std::size_t i = 0; i < breaks.size(); ++i) {
        breaks[i] = 2.0 * static_cast<double>(i) / (breaks.size() - 1);
    }

    NonUniformMeshX mesh(breaks, MCoordX(0));
    MDomainImpl<NonUniformMeshX> const dom(mesh, MCoordX(3));
    // EXPECT_EQ(mesh.lbound(), dom.lbound());
    // EXPECT_EQ(mesh.ubound(), dom.ubound());
    EXPECT_EQ(dom.rmin(), breaks.front());
    EXPECT_EQ(dom.rmax(), breaks.back());
}
