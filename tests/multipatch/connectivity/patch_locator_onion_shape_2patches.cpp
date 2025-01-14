// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_uniform.hpp"
#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
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
    using PatchLocator = OnionPatchLocator<
            MultipatchType<IdxRangeOnPatch, Patch1, Patch2>,
            CircularToCartesian<R, Theta, X, Y>,
            CartesianToCircular<X, Y, R, Theta>>;

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
        Patch1::Coord1 const r1_min(0.2);
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

template <class PatchLocator>
void test_operator_assignement(
        PatchLocator const& patch_locator,
        Kokkos::View<Coord<X, Y>*> coords,
        Kokkos::View<int*> patches)
{
    std::size_t failed_attempt = 0;
    Kokkos::parallel_reduce(
            coords.extent(0),
            KOKKOS_LAMBDA(int const i, std::size_t& attempt) {
                int patch_idx = patch_locator(coords(i));
                if (patch_idx != patches(i)) {
                    attempt++;
                }
            },
            failed_attempt);
    EXPECT_EQ(failed_attempt, 0);
}

} // namespace


TEST_F(OnionPatchLocator2PatchesTest, CheckPatchOrderingTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CzarnyToCartesian<R, Theta, X, Y> to_physical_mapping(0.3, 1.4);
    CartesianToCzarny<X, Y, R, Theta> to_logical_mapping(0.3, 1.4);

    EXPECT_NO_THROW((OnionPatchLocator(all_idx_ranges, to_physical_mapping, to_logical_mapping)));
}

TEST_F(OnionPatchLocator2PatchesTest, CheckPatchOrderingDeathTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch2, Patch1> all_idx_ranges(idx_range_2, idx_range_1);

    CircularToCartesian<R, Theta, X, Y> to_physical_mapping;
    CartesianToCircular<X, Y, R, Theta> to_logical_mapping;

    EXPECT_THROW(
            (OnionPatchLocator<
                    MultipatchType<IdxRangeOnPatch, Patch2, Patch1>,
                    CircularToCartesian<R, Theta, X, Y>,
                    CartesianToCircular<X, Y, R, Theta>,
                    Kokkos::DefaultExecutionSpace>(
                    all_idx_ranges,
                    to_physical_mapping,
                    to_logical_mapping)),
            std::invalid_argument);
}

TEST_F(OnionPatchLocator2PatchesTest, DeviceCircularOnionPatchLocator2PatchesTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CircularToCartesian<R, Theta, X, Y> to_physical_mapping;
    CartesianToCircular<X, Y, R, Theta> to_logical_mapping;

    OnionPatchLocator<
            MultipatchType<IdxRangeOnPatch, Patch1, Patch2>,
            CircularToCartesian<R, Theta, X, Y>,
            CartesianToCircular<X, Y, R, Theta>,
            Kokkos::DefaultExecutionSpace>
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);

    int constexpr n_elements = 7;
    Kokkos::View<Coord<X, Y>*> coords("coords 1", n_elements);
    Kokkos::View<int*> patches("patches 1", n_elements);

    Kokkos::View<Coord<X, Y>*, Kokkos::DefaultHostExecutionSpace>
            coords_host("coords_host 1", n_elements);
    Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace>
            patches_host("patches_host 1", n_elements);

    coords_host(0) = PhysicalCoordXY(0.15, .03);
    coords_host(1) = PhysicalCoordXY(0.2, 0);
    coords_host(2) = PhysicalCoordXY(0.25, .03);
    coords_host(3) = PhysicalCoordXY(1, 1);
    coords_host(4) = PhysicalCoordXY(1, 0);
    coords_host(5) = PhysicalCoordXY(1.5, 0);
    coords_host(6) = PhysicalCoordXY(-2.1, 0);

    patches_host(0) = PatchLocator::outside_rmin_domain;
    patches_host(1) = 0;
    patches_host(2) = 0;
    patches_host(3) = 1;
    patches_host(4) = 1;
    patches_host(5) = 1;
    patches_host(6) = PatchLocator::outside_rmax_domain;

    Kokkos::deep_copy(coords, coords_host);
    Kokkos::deep_copy(patches, patches_host);

    test_operator_assignement(patch_locator, coords, patches);
}


TEST_F(OnionPatchLocator2PatchesTest, DeviceCzarnyOnionPatchLocator2PatchesTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CzarnyToCartesian<R, Theta, X, Y> to_physical_mapping(0.3, 1.4);
    CartesianToCzarny<X, Y, R, Theta> to_logical_mapping(0.3, 1.4);

    OnionPatchLocator<
            MultipatchType<IdxRangeOnPatch, Patch1, Patch2>,
            CzarnyToCartesian<R, Theta, X, Y>,
            CartesianToCzarny<X, Y, R, Theta>,
            Kokkos::DefaultExecutionSpace>
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);

    int constexpr n_elements = 7;
    Kokkos::View<Coord<X, Y>*> coords("coords 1", n_elements);
    Kokkos::View<int*> patches("patches 1", n_elements);

    Kokkos::View<Coord<X, Y>*, Kokkos::DefaultHostExecutionSpace>
            coords_host("coords_host 1", n_elements);
    Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace>
            patches_host("patches_host 1", n_elements);

    coords_host(0) = PhysicalCoordXY(to_physical_mapping(Coord<R, Theta>(0, 0)));
    coords_host(1) = PhysicalCoordXY(to_physical_mapping(Coord<R, Theta>(0.2, 0)));
    coords_host(2) = PhysicalCoordXY(0.25, .03);
    coords_host(3) = PhysicalCoordXY(1, 1);
    coords_host(4) = PhysicalCoordXY(to_physical_mapping(Coord<R, Theta>(1, 0)));
    coords_host(5) = PhysicalCoordXY(to_physical_mapping(Coord<R, Theta>(1.5, 0)));
    coords_host(6) = PhysicalCoordXY(-2.1, 0);

    patches_host(0) = PatchLocator::outside_rmin_domain;
    patches_host(1) = 0;
    patches_host(2) = 0;
    patches_host(3) = 1;
    patches_host(4) = 1;
    patches_host(5) = 1;
    patches_host(6) = PatchLocator::outside_rmax_domain;

    Kokkos::deep_copy(coords, coords_host);
    Kokkos::deep_copy(patches, patches_host);

    test_operator_assignement(patch_locator, coords, patches);
}


TEST_F(OnionPatchLocator2PatchesTest, HostCircularOnionPatchLocator2PatchesTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CircularToCartesian<R, Theta, X, Y> to_physical_mapping;
    CartesianToCircular<X, Y, R, Theta> to_logical_mapping;
    using MultipatchIdxRanges = MultipatchType<IdxRangeOnPatch, Patch1, Patch2>;
    using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
    using PhysicalToLogicalMapping = CartesianToCircular<X, Y, R, Theta>;

    OnionPatchLocator<
            MultipatchIdxRanges,
            LogicalToPhysicalMapping,
            PhysicalToLogicalMapping,
            Kokkos::DefaultHostExecutionSpace>
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);

    PhysicalCoordXY coord(0.15, .03);
    EXPECT_EQ(patch_locator(coord), PatchLocator::outside_rmin_domain);

    coord = PhysicalCoordXY(0.25, .03);
    EXPECT_EQ(patch_locator(coord), 0);

    coord = PhysicalCoordXY(0.2, 0);
    EXPECT_EQ(patch_locator(coord), 0);

    coord = PhysicalCoordXY(1, 1);
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(1, 0);
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(1.5, 0);
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(-2.1, 0);
    EXPECT_EQ(patch_locator(coord), PatchLocator::outside_rmax_domain);
}


TEST_F(OnionPatchLocator2PatchesTest, HostCzarnyOnionPatchLocator2PatchesTest)
{
    Patch1::IdxRange12 idx_range_1(idx_range_r1, idx_range_theta1);
    Patch2::IdxRange12 idx_range_2(idx_range_r2, idx_range_theta2);
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(idx_range_1, idx_range_2);

    CzarnyToCartesian<R, Theta, X, Y> to_physical_mapping(0.3, 1.4);
    CartesianToCzarny<X, Y, R, Theta> to_logical_mapping(0.3, 1.4);

    using MultipatchIdxRanges = MultipatchType<IdxRangeOnPatch, Patch1, Patch2>;
    using LogicalToPhysicalMapping = CzarnyToCartesian<R, Theta, X, Y>;
    using PhysicalToLogicalMapping = CartesianToCzarny<X, Y, R, Theta>;
    OnionPatchLocator<
            MultipatchIdxRanges,
            LogicalToPhysicalMapping,
            PhysicalToLogicalMapping,
            Kokkos::DefaultHostExecutionSpace>
            patch_locator(all_idx_ranges, to_physical_mapping, to_logical_mapping);

    PhysicalCoordXY coord(-0.05, .03);
    EXPECT_EQ(patch_locator(coord), PatchLocator::outside_rmin_domain);

    coord = PhysicalCoordXY(to_physical_mapping(Coord<R, Theta>(0.2, 0)));
    EXPECT_EQ(patch_locator(coord), 0);

    coord = PhysicalCoordXY(0.25, .03);
    EXPECT_EQ(patch_locator(coord), 0);

    coord = PhysicalCoordXY(1, 1);
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(to_physical_mapping(Coord<R, Theta>(1, 0)));
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(to_physical_mapping(Coord<R, Theta>(1.5, 0)));
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(-2.1, 0);
    EXPECT_EQ(patch_locator(coord), PatchLocator::outside_rmax_domain);
}
