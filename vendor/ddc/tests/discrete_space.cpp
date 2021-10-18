// SPDX-License-Identifier: MIT
#include <array>
#include <iosfwd>
#include <memory>

#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/NonUniformDiscretization>
#include <ddc/UniformDiscretization>
#include <ddc/detail/discrete_space.hpp>

#include <gtest/gtest.h>

namespace {

class DimX;
class DimY;

using CoordX = Coordinate<DimX>;
using RCoordY = Coordinate<DimY>;
using DDimX = UniformDiscretization<DimX>;
using DDimY = NonUniformDiscretization<DimY>;

using IndexX = DiscreteCoordinate<DDimX>;

using DDimXY = detail::DiscreteSpace<DDimX, DDimY>;
using MCoordXY = DiscreteCoordinate<DDimX, DDimY>;
using RCoordXY = Coordinate<DimX, DimY>;

constexpr DDimX ddim_x = DDimX(CoordX(2.), CoordX(0.1));
constexpr std::array<RCoordY, 4> points_y {RCoordY(-1.), RCoordY(0.), RCoordY(2.), RCoordY(4.)};
static DDimY const ddim_y = DDimY(points_y);
constexpr DDimXY ddim_x_y {ddim_x, ddim_y};

} // namespace

TEST(DiscreteSpaceTest, constructor)
{
    EXPECT_EQ(DDimXY::rank(), DDimX::rank() + DDimY::rank());
    EXPECT_EQ(ddim_x_y.to_real(MCoordXY(0, 0)), RCoordXY(2., -1.));
}

TEST(DiscreteSpaceTest, accessor)
{
    EXPECT_EQ(ddim_x_y.get<DDimX>(), ddim_x);
}

TEST(DiscreteSpaceTest, submesh)
{
    auto&& selection = select<DDimY>(ddim_x_y);
    EXPECT_EQ(1, selection.rank());
    EXPECT_EQ(RCoordY(0.), selection.to_real(DiscreteCoordinate<DDimY>(1)));
}

TEST(DiscreteSpaceTest, conversion)
{
    constexpr static DDimX ddim_x(CoordX(2.), CoordX(0.1));
    constexpr detail::DiscreteSpace product_mesh_x(ddim_x);
    DDimX const& ddim_x_ref = get<DDimX>(product_mesh_x);
    double step = ddim_x_ref.step();
    EXPECT_EQ(0.1, step);
}

TEST(DiscreteSpaceTest, to_real)
{
    for (std::size_t ix = 0; ix < 5; ++ix) {
        for (std::size_t ivx = 0; ivx < points_y.size(); ++ivx) {
            EXPECT_EQ(
                    RCoordXY(ddim_x.to_real(IndexX(ix)), points_y[ivx]),
                    ddim_x_y.to_real(MCoordXY(ix, ivx)));
        }
    }
}
