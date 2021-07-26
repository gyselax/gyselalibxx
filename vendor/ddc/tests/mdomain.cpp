#include <iosfwd>
#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "gtest/gtest_pred_impl.h"

#include "mcoord.h"
#include "mdomain.h"
#include "rcoord.h"
#include "uniform_mesh.h"

class DimX;

using MeshX = UniformMesh<DimX>;
using MCoordX = MCoord<MeshX>;

using RCoordX = RCoord<DimX>;

class MDomainXTest : public ::testing::Test
{
protected:
    std::size_t npoints = 101;
    MeshX mesh_x = MeshX(0., 1., npoints);
    MDomain<MeshX> const dom = MDomain(mesh_x, npoints - 1);
};

TEST_F(MDomainXTest, Constructor)
{
    constexpr RCoordX origin(1);
    constexpr RCoordX unit_vec(3);
    constexpr MCoordX lbound(0);
    constexpr MCoordX ubound(10);
    constexpr RCoordX rmin = origin;
    constexpr RCoordX rmax = rmin + unit_vec * (double)(ubound - lbound);
    constexpr static MeshX mesh(origin, unit_vec);
    constexpr MDomain dom_a(mesh, ubound);
    constexpr MDomain dom_b(mesh, lbound, ubound);
    constexpr MDomain dom_d(dom_a);
    EXPECT_EQ(dom_a, dom_b);
    EXPECT_EQ(dom_a, dom_d);
    EXPECT_EQ(dom_a.rmin(), rmin);
    EXPECT_EQ(dom_a.rmax(), rmax);
    EXPECT_EQ(dom_a.lbound(), lbound);
    EXPECT_EQ(dom_a.ubound(), ubound);
}

TEST_F(MDomainXTest, subdomain)
{
    EXPECT_EQ(dom.lbound(), 0ul);
    EXPECT_EQ(dom.ubound(), 100ul);
    EXPECT_EQ(dom.size(), npoints);

    MDomain const subdomain1 = dom.subdomain(10ul, 91ul);
    EXPECT_EQ(subdomain1.lbound(), 10ul);
    EXPECT_EQ(subdomain1.ubound(), 100ul);
    EXPECT_EQ(subdomain1.size(), 91ul);

    MDomain const subdomain2 = subdomain1.subdomain(10ul, 2ul);
    EXPECT_EQ(subdomain2.lbound(), 20ul);
    EXPECT_EQ(subdomain2.ubound(), 21ul);
    EXPECT_EQ(subdomain2.size(), 2ul);
}

TEST_F(MDomainXTest, ubound)
{
    EXPECT_EQ(dom.ubound().get<MeshX>(), 100ul);
    EXPECT_EQ(dom.ubound(), 100ul);
}

TEST_F(MDomainXTest, rmax)
{
    EXPECT_EQ(dom.mesh().to_real(get<MeshX>(dom.ubound())).get<DimX>(), 1.);
    EXPECT_EQ(dom.rmax().get<DimX>(), 1.);
    EXPECT_EQ(dom.rmax(), 1.);
}

TEST_F(MDomainXTest, RangeFor)
{
    std::size_t ii = 0;
    for (auto&& x : dom) {
        ASSERT_LE(0, x);
        EXPECT_EQ(x, ii);
        ASSERT_LT(x, npoints);
        ++ii;
    }
}
