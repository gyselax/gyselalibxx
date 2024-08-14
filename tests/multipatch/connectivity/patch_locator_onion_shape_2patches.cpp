// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_uniform.hpp"
#include "multipatch_type.hpp"
#include "onion_patch_locator.hpp"
#include "patch.hpp"
#include "physical_geometry.hpp"
#include "types.hpp"


using namespace onion_shape_uniform_2d_2patches;
using namespace physical_geometry;


namespace {
class OnionPatchLocator2PatchesTest : public ::testing::Test
{
private:
    static constexpr Patch1::IdxStep1 r1_size = Patch1::IdxStep1(10);
    static constexpr Patch1::IdxStep2 theta1_size = Patch1::IdxStep2(10);

    static constexpr Patch2::IdxStep1 r2_size = Patch2::IdxStep1(10);
    static constexpr Patch2::IdxStep2 theta2_size = Patch2::IdxStep2(10);

public:
    static constexpr int outside_domain = IPatchLocator<X, Y>::outside_domain;

    Patch1::IdxRange1 const idx_range_r1;
    Patch1::IdxRange2 const idx_range_theta1;

    Patch2::IdxRange1 const idx_range_r2;
    Patch2::IdxRange2 const idx_range_theta2;

public:
    OnionPatchLocator2PatchesTest()
        : idx_range_r1(Patch1::Idx1(0), r1_size)
        , idx_range_theta1(Patch1::Idx2(0), theta1_size)
        , idx_range_r2(Patch2::Idx1(0), r2_size)
        , idx_range_theta2(Patch2::Idx2(0), theta2_size) {};

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        Patch1::Coord1 const r1_min(0.0);
        Patch1::Coord1 const r1_max(1.0);

        Patch1::Coord2 const theta1_min(0.0);
        Patch1::Coord2 const theta1_max(2 * M_PI);

        ddc::init_discrete_space<GridR<1>>(GridR<1>::init(r1_min, r1_max, r1_size));
        ddc::init_discrete_space<GridTheta<1>>(
                GridTheta<1>::init(theta1_min, theta1_max, theta1_size));

        // Patch 2
        Patch2::Coord1 const r2_min(1.0);
        Patch2::Coord1 const r2_max(1.5);

        Patch2::Coord2 const theta2_min(0.0);
        Patch2::Coord2 const theta2_max(2 * M_PI);

        ddc::init_discrete_space<GridR<2>>(GridR<2>::init(r2_min, r2_max, r2_size));
        ddc::init_discrete_space<GridTheta<2>>(
                GridTheta<2>::init(theta2_min, theta2_max, theta2_size));
    }
};

} // namespace


TEST_F(OnionPatchLocator2PatchesTest, CheckPatchOrderingTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CzarnyToCartesian<X, Y, R, Theta> mapping(0.3, 1.4);

    EXPECT_NO_THROW(OnionPatchLocator(all_idx_ranges, mapping));
}

TEST_F(OnionPatchLocator2PatchesTest, CheckPatchOrderingDeathTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch2, Patch1> all_idx_ranges(idx_range_2, idx_range_1);

    CircularToCartesian<X, Y, R, Theta> mapping;

    EXPECT_THROW(OnionPatchLocator(all_idx_ranges, mapping), std::invalid_argument);
}



TEST_F(OnionPatchLocator2PatchesTest, CircularOnionPatchLocator2PatchesTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CircularToCartesian<X, Y, R, Theta> mapping;

    OnionPatchLocator patch_locator(all_idx_ranges, mapping);

    PhysicalCoordXY coord(0.25, .03);
    int patch = patch_locator(coord);
    EXPECT_EQ(patch, 0);

    coord = PhysicalCoordXY(1, 1);
    patch = patch_locator(coord);
    EXPECT_EQ(patch, 1);

    coord = PhysicalCoordXY(-2.1, 0);
    patch = patch_locator(coord);
    EXPECT_EQ(patch, outside_domain);
}


TEST_F(OnionPatchLocator2PatchesTest, CzarnyOnionPatchLocator2PatchesTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CzarnyToCartesian<X, Y, R, Theta> mapping(0.3, 1.4);

    OnionPatchLocator patch_locator(all_idx_ranges, mapping);

    PhysicalCoordXY coord(0.25, .03);
    int patch = patch_locator(coord);
    EXPECT_EQ(patch, 0);

    coord = PhysicalCoordXY(1, 1);
    patch = patch_locator(coord);
    EXPECT_EQ(patch, 1);

    coord = PhysicalCoordXY(-2.1, 0);
    patch = patch_locator(coord);
    EXPECT_EQ(patch, outside_domain);
}