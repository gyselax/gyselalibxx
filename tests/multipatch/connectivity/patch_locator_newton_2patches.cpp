// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_uniform.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "l_norm_tools.hpp"
#include "multipatch_type.hpp"
#include "patch.hpp"
#include "patch_locator_newton.hpp"
#include "physical_geometry.hpp"
#include "types.hpp"

using namespace onion_shape_uniform_2d_2patches;
using namespace physical_geometry;


namespace {

template <class P>
using CzarnyOnPatch = CzarnyToCartesian<R, Theta, X, Y>;

template <class P>
using CircularOnPatch = CircularToCartesian<R, Theta, X, Y>;
;

using MPCzarnyMapping = MultipatchType<CzarnyOnPatch, Patch1, Patch2>;
using MPCircularMapping = MultipatchType<CircularOnPatch, Patch1, Patch2>;

class PatchLocatorNewton2PatchesTest : public ::testing::Test
{
private:
    static constexpr Patch1::IdxStep1 r1_size = Patch1::IdxStep1(10);
    static constexpr Patch1::IdxStep2 theta1_size = Patch1::IdxStep2(10);

    static constexpr Patch2::IdxStep1 r2_size = Patch2::IdxStep1(10);
    static constexpr Patch2::IdxStep2 theta2_size = Patch2::IdxStep2(10);

public:
    static constexpr Patch1::Coord1 r1_min = Patch1::Coord1(0.2);
    static constexpr Patch1::Coord1 r1_max = Patch1::Coord1(1.0);

    static constexpr Patch1::Coord2 theta1_min = Patch1::Coord2(0.0);
    static constexpr Patch1::Coord2 theta1_max = Patch1::Coord2(2 * M_PI);

    static constexpr Patch2::Coord1 r2_min = Patch2::Coord1(1.0);
    static constexpr Patch2::Coord1 r2_max = Patch2::Coord1(1.5);

    static constexpr Patch2::Coord2 theta2_min = Patch2::Coord2(0.0);
    static constexpr Patch2::Coord2 theta2_max = Patch2::Coord2(2 * M_PI);

    Patch1::IdxRange1 const idx_range_r1;
    Patch1::IdxRange2 const idx_range_theta1;

    Patch2::IdxRange1 const idx_range_r2;
    Patch2::IdxRange2 const idx_range_theta2;

public:
    PatchLocatorNewton2PatchesTest()
        : idx_range_r1(Patch1::Idx1(0), r1_size)
        , idx_range_theta1(Patch1::Idx2(0), theta1_size)
        , idx_range_r2(Patch2::Idx1(0), r2_size)
        , idx_range_theta2(Patch2::Idx2(0), theta2_size) {};

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        ddc::init_discrete_space<BSplinesR<1>>(r1_min, r1_max, r1_size - 1);
        ddc::init_discrete_space<BSplinesTheta<1>>(theta1_min, theta1_max, theta1_size - 1);

        ddc::init_discrete_space<GridR<1>>(GridR<1>::init(r1_min, r1_max, r1_size));
        ddc::init_discrete_space<GridTheta<1>>(
                GridTheta<1>::init(theta1_min, theta1_max, theta1_size));

        // Patch 2
        ddc::init_discrete_space<BSplinesR<2>>(r2_min, r2_max, r2_size - 1);
        ddc::init_discrete_space<BSplinesTheta<2>>(theta2_min, theta2_max, theta2_size - 1);

        ddc::init_discrete_space<GridR<2>>(GridR<2>::init(r2_min, r2_max, r2_size));
        ddc::init_discrete_space<GridTheta<2>>(
                GridTheta<2>::init(theta2_min, theta2_max, theta2_size));
    }
};

template <class PatchLocator>
void test_patch_index(
        PatchLocator const& patch_locator,
        Kokkos::View<Coord<X, Y>*> coords,
        Kokkos::View<int*> patches)
{
    std::size_t failed_attempt = 0;
    Kokkos::parallel_reduce(
            coords.extent(0),
            KOKKOS_LAMBDA(int const i, std::size_t& attempt) {
                int patch_idx = patch_locator(coords(i));
                std::cout << patches(i) << "  " << patch_idx << "  " << coords(i) << std::endl;
                if (patch_idx != patches(i)) {
                    attempt++;
                }
            },
            failed_attempt);
    EXPECT_EQ(failed_attempt, 0);
}

template <class PatchLocator, class Mapping>
void test_transformation_inversion(
        PatchLocator const& patch_locator,
        Kokkos::View<Coord<R, Theta>*> coords_log,
        Mapping const& mapping)
{
    std::size_t failed_attempt = 0;
    Kokkos::parallel_reduce(
            coords_log.extent(0),
            KOKKOS_LAMBDA(int const i, std::size_t& attempt) {
                Coord<X, Y> coord = mapping(coords_log(i));
                Coord<R, Theta> coord_log = coords_log(i) + Coord<R, Theta>(0.3, -0.2);
                int patch_idx = patch_locator(coord_log, coord);
                double err_r = ddc::get<R>(coord_log - coords_log(i));
                double delta_theta = Kokkos::
                        fmod(Kokkos::fabs(ddc::get<Theta>(coord_log - coords_log(i))), 2 * M_PI);
                double err_theta = Kokkos::fmin(delta_theta, 2 * M_PI - delta_theta);
                if ((patch_idx >= 0) && ((err_r > 1e-12) || (err_theta > 1e-11))) {
                    attempt++;
                }
            },
            failed_attempt);
    EXPECT_EQ(failed_attempt, 0);
}

template <class PatchLocator, class Mapping>
void test_transformation_inversion_with_init(
        PatchLocator const& patch_locator,
        Kokkos::View<Coord<R, Theta>*> coords_log,
        Mapping const& mapping)
{
    std::size_t failed_attempt = 0;
    Kokkos::parallel_reduce(
            coords_log.extent(0),
            KOKKOS_LAMBDA(int const i, std::size_t& attempt) {
                Coord<X, Y> coord = mapping(coords_log(i));
                Coord<R, Theta> coord_log = coords_log(i) + Coord<R, Theta>(0.3, -0.2);
                int patch_idx = patch_locator.template find<Patch1>(coord_log, coord);
                double err_r = ddc::get<R>(coord_log - coords_log(i));
                double delta_theta = Kokkos::
                        fmod(Kokkos::fabs(ddc::get<Theta>(coord_log - coords_log(i))), 2 * M_PI);
                double err_theta = Kokkos::fmin(delta_theta, 2 * M_PI - delta_theta);
                if ((patch_idx >= 0) && ((err_r > 1e-12) || (err_theta > 1e-11))) {
                    attempt++;
                }
            },
            failed_attempt);
    EXPECT_EQ(failed_attempt, 0);
}

} // namespace


TEST_F(PatchLocatorNewton2PatchesTest, DeviceCircularTest)
{
    using PatchLocator = PatchLocatorNewton<MPCircularMapping>;

    CircularOnPatch<Patch1> mapping_1;
    CircularOnPatch<Patch2> mapping_2;
    MPCircularMapping mappings(mapping_1, mapping_2);

    PatchLocator patch_locator(mappings);

    int constexpr n_elements = 7;
    Kokkos::View<Coord<X, Y>*> coords("coords", n_elements);
    Kokkos::View<int*> patches("patches", n_elements);

    Kokkos::View<Coord<X, Y>*, Kokkos::DefaultHostExecutionSpace>
            coords_host("coords_host", n_elements);
    Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> patches_host("patches_host", n_elements);

    coords_host(0) = PhysicalCoordXY(0.15, .03);
    coords_host(1) = PhysicalCoordXY(0.2, 0);
    coords_host(2) = PhysicalCoordXY(0.25, .03);
    coords_host(3) = PhysicalCoordXY(1, 1);
    coords_host(4) = PhysicalCoordXY(1, 0);
    coords_host(5) = PhysicalCoordXY(1.5, 0);
    coords_host(6) = PhysicalCoordXY(-2.1, 0);

    patches_host(0) = PatchLocator::outside_domain;
    patches_host(1) = 0;
    patches_host(2) = 0;
    patches_host(3) = 1;
    patches_host(4) = 0;
    patches_host(5) = 1;
    patches_host(6) = PatchLocator::outside_domain;

    Kokkos::deep_copy(coords, coords_host);
    Kokkos::deep_copy(patches, patches_host);

    test_patch_index(patch_locator, coords, patches);
}


TEST_F(PatchLocatorNewton2PatchesTest, DeviceCzarnyTest)
{
    using PatchLocator = PatchLocatorNewton<MPCzarnyMapping>;

    CzarnyOnPatch<Patch1> mapping_1(0.3, 1.4);
    CzarnyOnPatch<Patch2> mapping_2(0.3, 1.4);
    MPCzarnyMapping mappings(mapping_1, mapping_2);

    PatchLocator patch_locator(mappings);

    CartesianToCzarny<X, Y, R, Theta> to_logical_mapping(0.3, 1.4);

    int constexpr n_elements = 7;
    Kokkos::View<Coord<X, Y>*> coords("coords", n_elements);
    Kokkos::View<int*> patches("patches", n_elements);

    Kokkos::View<Coord<X, Y>*, Kokkos::DefaultHostExecutionSpace>
            coords_host("coords_host", n_elements);
    Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> patches_host("patches_host", n_elements);

    // Test outside domain (radius smaller than the minimum radius)
    coords_host(0) = PhysicalCoordXY(mapping_1(Coord<R, Theta>(0, 0)));
    // Test radial lower bound
    coords_host(1) = PhysicalCoordXY(mapping_1(Coord<R, Theta>(r1_min, theta1_min)));
    // Test in patch 1
    coords_host(2) = PhysicalCoordXY(0.25, .03);
    // Test in patch 2
    coords_host(3) = PhysicalCoordXY(1, 1);
    // Test patch interface
    coords_host(4) = PhysicalCoordXY(mapping_1(Coord<R, Theta>(r1_max, theta1_min)));
    // Test radial upper bound
    coords_host(5) = PhysicalCoordXY(mapping_1(Coord<R, Theta>(r2_max, theta2_min)));
    // Test outside domain (radius bigger than the maximum radius)
    coords_host(6) = PhysicalCoordXY(-2.1, 0);

    // Test outside domain (radius smaller than the minimum radius)
    patches_host(0) = PatchLocator::outside_domain;
    // Test radial lower bound
    if (ddc::select<R>(to_logical_mapping(coords_host(1))) - r1_min < 0) {
        patches_host(1) = PatchLocator::outside_domain;
    } else {
        patches_host(1) = 0;
    }
    // Test in patch 1
    patches_host(2) = 0;
    // Test in patch 2
    patches_host(3) = 1;
    // Test patch interface
    if (ddc::select<R>(to_logical_mapping(coords_host(4))) - r1_max <= 1e-12) {
        patches_host(4) = 0;
    } else {
        patches_host(4) = 1;
    }
    // Test radial upper bound
    if (ddc::select<R>(to_logical_mapping(coords_host(5))) - r2_max <= 1e-12) {
        patches_host(5) = 1;
    } else {
        patches_host(5) = PatchLocator::outside_domain;
    }
    // Test outside domain (radius bigger than the maximum radius)
    patches_host(6) = PatchLocator::outside_domain;

    Kokkos::deep_copy(coords, coords_host);
    Kokkos::deep_copy(patches, patches_host);

    test_patch_index(patch_locator, coords, patches);
}


TEST_F(PatchLocatorNewton2PatchesTest, HostCircularTest)
{
    using PatchLocator = PatchLocatorNewton<MPCircularMapping>;

    CircularOnPatch<Patch1> mapping_1;
    CircularOnPatch<Patch2> mapping_2;
    MPCircularMapping mappings(mapping_1, mapping_2);

    PatchLocator patch_locator(mappings);

    PhysicalCoordXY coord(0.15, .03);
    EXPECT_EQ(patch_locator(coord), PatchLocator::outside_domain);

    coord = PhysicalCoordXY(0.25, .03);
    EXPECT_EQ(patch_locator(coord), 0);

    coord = PhysicalCoordXY(0.2, 0);
    EXPECT_EQ(patch_locator(coord), 0);

    coord = PhysicalCoordXY(1, 1);
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(1, 0);
    EXPECT_EQ(patch_locator(coord), 0);

    coord = PhysicalCoordXY(1.5, 0);
    EXPECT_EQ(patch_locator(coord), 1);

    coord = PhysicalCoordXY(-2.1, 0);
    EXPECT_EQ(patch_locator(coord), PatchLocator::outside_domain);
}

TEST_F(PatchLocatorNewton2PatchesTest, DeviceCircularInversionTest)
{
    using PatchLocator = PatchLocatorNewton<MPCircularMapping>;

    CircularOnPatch<Patch1> mapping_1;
    CircularOnPatch<Patch2> mapping_2;
    MPCircularMapping mappings(mapping_1, mapping_2);

    PatchLocator patch_locator(mappings);

    int constexpr n_elements = 7;
    Kokkos::View<Coord<R, Theta>*> coords_log("coords_logical", n_elements);
    Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultHostExecutionSpace>
            coords_log_host("coords_logical_host", n_elements);

    coords_log_host(0) = Coord<R, Theta>(0.2, 0);
    coords_log_host(1) = Coord<R, Theta>(1.0, M_PI);
    coords_log_host(2) = Coord<R, Theta>(1.5, 2 * M_PI);
    coords_log_host(3) = Coord<R, Theta>(0, 0.5 * M_PI);
    coords_log_host(4) = Coord<R, Theta>(0.99, M_PI / 4);
    coords_log_host(5) = Coord<R, Theta>(1.0001, M_PI / 2);
    coords_log_host(6) = Coord<R, Theta>(1.49, M_PI * 3 / 2);

    Kokkos::deep_copy(coords_log, coords_log_host);

    test_transformation_inversion(patch_locator, coords_log, mapping_1);
}

TEST_F(PatchLocatorNewton2PatchesTest, DeviceCzarnyInversionTest)
{
    using PatchLocator = PatchLocatorNewton<MPCzarnyMapping>;

    CzarnyOnPatch<Patch1> mapping_1(0.3, 1.4);
    CzarnyOnPatch<Patch2> mapping_2(0.3, 1.4);
    MPCzarnyMapping mappings(mapping_1, mapping_2);

    PatchLocator patch_locator(mappings);

    int constexpr n_elements = 7;
    Kokkos::View<Coord<R, Theta>*> coords_log("coords_logical", n_elements);
    Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultHostExecutionSpace>
            coords_log_host("coords_logical_host", n_elements);

    coords_log_host(0) = Coord<R, Theta>(0.2, 0);
    coords_log_host(1) = Coord<R, Theta>(1.0, M_PI);
    coords_log_host(2) = Coord<R, Theta>(1.5, 2 * M_PI);
    coords_log_host(3) = Coord<R, Theta>(0, 0.5 * M_PI);
    coords_log_host(4) = Coord<R, Theta>(0.99, M_PI / 4);
    coords_log_host(5) = Coord<R, Theta>(1.0001, M_PI / 2);
    coords_log_host(6) = Coord<R, Theta>(1.49, M_PI * 3 / 2);

    Kokkos::deep_copy(coords_log, coords_log_host);

    test_transformation_inversion(patch_locator, coords_log, mapping_1);
}

TEST_F(PatchLocatorNewton2PatchesTest, DeviceCzarnyInversionWithInitTest)
{
    using PatchLocator = PatchLocatorNewton<MPCzarnyMapping>;

    CzarnyOnPatch<Patch1> mapping_1(0.3, 1.4);
    CzarnyOnPatch<Patch2> mapping_2(0.3, 1.4);
    MPCzarnyMapping mappings(mapping_1, mapping_2);

    PatchLocator patch_locator(mappings);

    int constexpr n_elements = 7;
    Kokkos::View<Coord<R, Theta>*> coords_log("coords_logical", n_elements);
    Kokkos::View<Coord<R, Theta>*, Kokkos::DefaultHostExecutionSpace>
            coords_log_host("coords_logical_host", n_elements);

    coords_log_host(0) = Coord<R, Theta>(0.2, 0);
    coords_log_host(1) = Coord<R, Theta>(1.0, M_PI);
    coords_log_host(2) = Coord<R, Theta>(1.5, 2 * M_PI);
    coords_log_host(3) = Coord<R, Theta>(0, 0.5 * M_PI);
    coords_log_host(4) = Coord<R, Theta>(0.99, M_PI / 4);
    coords_log_host(5) = Coord<R, Theta>(1.0001, M_PI / 2);
    coords_log_host(6) = Coord<R, Theta>(1.49, M_PI * 3 / 2);

    Kokkos::deep_copy(coords_log, coords_log_host);

    test_transformation_inversion_with_init(patch_locator, coords_log, mapping_1);
}