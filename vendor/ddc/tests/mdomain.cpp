// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>
#include <utility>

#include <ddc/MCoord>
#include <ddc/MDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>

#include <gtest/gtest.h>

class DimX;

using MeshX = UniformMesh<DimX>;
using MCoordX = MCoord<MeshX>;
using MLengthX = MLength<MeshX>;

using RCoordX = RCoord<DimX>;

class MDomainXTest : public ::testing::Test
{
protected:
    MeshX mesh_x = MeshX(RCoordX(0.), RCoordX(0.01));
    MCoordX lbound_x = 50;
    MLengthX size_x = 51;
    MDomain<MeshX> const dom = MDomain(mesh_x, lbound_x, size_x);
    RCoordX min_x = 0.5;
    RCoordX max_x = 1.;
    MCoordX ubound_x = 100;
};

TEST_F(MDomainXTest, Constructor)
{
    EXPECT_EQ(dom.mesh(), mesh_x);
    EXPECT_EQ(dom.size(), size_x);
    EXPECT_EQ(dom.empty(), false);
    EXPECT_EQ(dom[0], lbound_x);
    EXPECT_EQ(dom.front(), lbound_x);
    EXPECT_EQ(dom.back(), ubound_x);
    EXPECT_EQ(dom.rmin(), min_x);
    EXPECT_EQ(dom.rmax(), max_x);
}

TEST_F(MDomainXTest, Empty)
{
    MLengthX size_x(0);
    MDomain const empty_domain(mesh_x, lbound_x, size_x);
    EXPECT_EQ(empty_domain.mesh(), mesh_x);
    EXPECT_EQ(empty_domain.size(), size_x);
    EXPECT_EQ(empty_domain.empty(), true);
    EXPECT_EQ(empty_domain[0], lbound_x);
}

TEST_F(MDomainXTest, RangeFor)
{
    MCoordX ii = lbound_x;
    for (MCoordX ix : dom) {
        ASSERT_LE(lbound_x, ix);
        EXPECT_EQ(ix, ii);
        ASSERT_LE(ix, ubound_x);
        ++ii.get<MeshX>();
    }
}
