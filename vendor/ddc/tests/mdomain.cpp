#include <gtest/gtest.h>

#include "mdomain.h"

class DimX;
class DimVx;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<DimX>;
using MeshX = UniformMesh<DimX>;
using MDomainX = UniformMDomain<DimX>;
using MeshXVx = UniformMesh<DimX, DimVx>;
using RCoordXVx = RCoord<DimX, DimVx>;
using MCoordXVx = MCoord<DimX, DimVx>;
using MDomainXVx = UniformMDomain<DimX, DimVx>;

TEST(MDomainXTest, Constructor)
{
    constexpr RCoordX origin(1);
    constexpr RCoordX unit_vec(3);
    constexpr MCoordX lbound(0);
    constexpr MCoordX ubound(10);
    constexpr RCoordX rmin = origin;
    constexpr RCoordX rmax = rmin + unit_vec * (ubound - lbound);
    constexpr MeshX mesh2d(origin, unit_vec);
    constexpr MDomainX dom2d_a(mesh2d, ubound);
    constexpr MDomainX dom2d_b(mesh2d, lbound, ubound);
    constexpr MDomainX dom2d_c(rmin, rmax, lbound, ubound);
    constexpr MDomainX dom2d_d(dom2d_a);
    EXPECT_EQ(dom2d_a, dom2d_b);
    EXPECT_EQ(dom2d_a, dom2d_c);
    EXPECT_EQ(dom2d_a, dom2d_d);
    EXPECT_EQ(dom2d_a.rmin(), rmin);
    EXPECT_EQ(dom2d_a.rmax(), rmax);
    EXPECT_EQ(dom2d_a.lbound(), lbound);
    EXPECT_EQ(dom2d_a.ubound(), ubound);
}

TEST(MDomainXVxTest, Constructor)
{
    constexpr RCoordXVx origin(1, 2);
    constexpr RCoordXVx unit_vec(3, 4);
    constexpr MCoordXVx lbound(0, 0);
    constexpr MCoordXVx ubound(10, 20);
    constexpr RCoordXVx rmin = origin;
    constexpr RCoordXVx rmax = rmin + unit_vec * (ubound - lbound);
    constexpr MeshXVx mesh2d(origin, unit_vec);
    constexpr MDomainXVx dom2d_a(mesh2d, ubound);
    constexpr MDomainXVx dom2d_b(mesh2d, lbound, ubound);
    constexpr MDomainXVx dom2d_c(rmin, rmax, lbound, ubound);
    constexpr MDomainXVx dom2d_d(dom2d_a);
    EXPECT_EQ(dom2d_a, dom2d_b);
    EXPECT_EQ(dom2d_a, dom2d_c);
    EXPECT_EQ(dom2d_a, dom2d_d);
    EXPECT_EQ(dom2d_a.rmin(), rmin);
    EXPECT_EQ(dom2d_a.rmax(), rmax);
    EXPECT_EQ(dom2d_a.lbound(), lbound);
    EXPECT_EQ(dom2d_a.ubound(), ubound);
}

TEST(MDomainXVxTest, ubound)
{
    MDomainXVx const dom2d(RCoordXVx(0, 0), RCoordXVx(2, 2), MCoordXVx(0, 0), MCoordXVx(100, 200));
    EXPECT_EQ(dom2d.ubound().get<DimX>(), 100ul);
    EXPECT_EQ(dom2d.ubound<DimX>(), 100ul);

    EXPECT_EQ(dom2d.ubound().get<DimVx>(), 200ul);
    EXPECT_EQ(dom2d.ubound<DimVx>(), 200ul);
}

TEST(MDomainXVxTest, rmax)
{
    MDomainXVx const
            dom2d(RCoordXVx(0., 200.),
                  RCoordXVx(200., 400.),
                  MCoordXVx(0ul, 0ul),
                  MCoordXVx(100ul, 200ul));
    EXPECT_EQ(dom2d.mesh().to_real(dom2d.ubound()).get<DimX>(), 200.);
    EXPECT_EQ(dom2d.rmax().get<DimX>(), 200.);
    EXPECT_EQ(dom2d.rmax<DimX>(), 200.);

    EXPECT_EQ(dom2d.mesh().to_real(dom2d.ubound()).get<DimVx>(), 400.);
    EXPECT_EQ(dom2d.rmax().get<DimVx>(), 400.);
    EXPECT_EQ(dom2d.rmax<DimVx>(), 400.);
}

TEST(MDomainXTest, RangeFor)
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
