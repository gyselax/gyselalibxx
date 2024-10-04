// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_non_uniform.hpp"
#include "constant_extrapolation_rules_onion.hpp"
#include "ddc_helper.hpp"
#include "mesh_builder.hpp"
#include "multipatch_type.hpp"
#include "onion_patch_locator.hpp"
#include "physical_geometry.hpp"
#include "types.hpp"

using namespace onion_shape_non_uniform_2d_2patches;
using namespace physical_geometry;


namespace {
using HostExecSpace = Kokkos::DefaultHostExecutionSpace;
using DeviceExecSpace = Kokkos::DefaultExecutionSpace;

template <int PatchIdx>
using SplineInterpPointsR = ddc::GrevilleInterpolationPoints<
        BSplinesR<PatchIdx>,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE>;
template <int PatchIdx>
using SplineInterpPointsTheta = ddc::GrevilleInterpolationPoints<
        BSplinesTheta<PatchIdx>,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;


template <int PatchIdx, class ExecSpace>
using SplineRThetaBuilder = ddc::SplineBuilder2D<
        ExecSpace,
        typename ExecSpace::memory_space,
        BSplinesR<PatchIdx>,
        BSplinesTheta<PatchIdx>,
        GridR<PatchIdx>,
        GridTheta<PatchIdx>,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridR<PatchIdx>,
        GridTheta<PatchIdx>>;


using Mapping = CircularToCartesian<X, Y, R, Theta>;
using MultipatchIdxRange = MultipatchType<IdxRangeOnPatch, Patch1, Patch2>;

template <class ExecSpace>
using PatchLocator = OnionPatchLocator<MultipatchIdxRange, Mapping, ExecSpace>;


static constexpr Coord<R> r_interface = Coord<R>(1.0);

class MultipatchConstantExtrapolationRuleOnionTest : public ::testing::Test
{
private:
    static constexpr Patch1::IdxStep1 r1_ncells = Patch1::IdxStep1(16);
    static constexpr Patch1::IdxStep2 theta1_ncells = Patch1::IdxStep2(10);

    static constexpr Patch2::IdxStep1 r2_ncells = Patch2::IdxStep1(8);
    static constexpr Patch2::IdxStep2 theta2_ncells = Patch2::IdxStep2(12);

protected:
    // Coordinates delimiting the patches
    static constexpr typename Patch1::Coord1 r1_min = Patch1::Coord1(0.25);
    static constexpr typename Patch1::Coord1 r1_max = r_interface;

    static constexpr typename Patch1::Coord2 theta1_min = Patch1::Coord2(0.0);
    static constexpr typename Patch1::Coord2 theta1_max = Patch1::Coord2(2 * M_PI);

    static constexpr typename Patch2::Coord1 r2_min = r_interface;
    static constexpr typename Patch2::Coord1 r2_max = Patch2::Coord1(2.0);

    static constexpr Patch2::Coord2 theta2_min = Patch2::Coord2(0.0);
    static constexpr Patch2::Coord2 theta2_max = Patch2::Coord2(2 * M_PI);

    // Index ranges
    Patch1::IdxRange12 const idx_range_rtheta1;
    Patch1::IdxRangeBS12 const spline_idx_range_rtheta1;
    Patch2::IdxRange12 const idx_range_rtheta2;
    Patch2::IdxRangeBS12 const spline_idx_range_rtheta2;

    // Spline representations
    DFieldMem<IdxRange<BSplinesR<1>, BSplinesTheta<1>>> function_1_coef_alloc;
    DFieldMem<IdxRange<BSplinesR<2>, BSplinesTheta<2>>> function_2_coef_alloc;

    MultipatchType<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const splines;

    MultipatchIdxRange const all_idx_ranges;


public:
    MultipatchConstantExtrapolationRuleOnionTest()
        : idx_range_rtheta1(
                SplineInterpPointsR<1>::get_domain<GridR<1>>(),
                SplineInterpPointsTheta<1>::get_domain<GridTheta<1>>())
        , spline_idx_range_rtheta1(
                  ddc::discrete_space<BSplinesR<1>>().full_domain(),
                  ddc::discrete_space<BSplinesTheta<1>>().full_domain())
        , idx_range_rtheta2(
                  SplineInterpPointsR<2>::get_domain<GridR<2>>(),
                  SplineInterpPointsTheta<2>::get_domain<GridTheta<2>>())
        , spline_idx_range_rtheta2(
                  ddc::discrete_space<BSplinesR<2>>().full_domain(),
                  ddc::discrete_space<BSplinesTheta<2>>().full_domain())
        , function_1_coef_alloc(spline_idx_range_rtheta1)
        , function_2_coef_alloc(spline_idx_range_rtheta2)
        , splines(set_spline_1<DeviceExecSpace>(get_field(function_1_coef_alloc)),
                  set_spline_2<DeviceExecSpace>(get_field(function_2_coef_alloc)))
        , all_idx_ranges(idx_range_rtheta1, idx_range_rtheta2)
    {
    }

    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        ddc::init_discrete_space<BSplinesR<1>>(
                build_uniform_break_points(r1_min, r1_max, r1_ncells));
        ddc::init_discrete_space<BSplinesTheta<1>>(
                build_uniform_break_points(theta1_min, theta1_max, theta1_ncells));

        ddc::init_discrete_space<GridR<1>>(SplineInterpPointsR<1>::get_sampling<GridR<1>>());
        ddc::init_discrete_space<GridTheta<1>>(
                SplineInterpPointsTheta<1>::get_sampling<GridTheta<1>>());

        // Patch 2
        ddc::init_discrete_space<BSplinesR<2>>(
                build_uniform_break_points(r2_min, r2_max, r2_ncells));
        ddc::init_discrete_space<BSplinesTheta<2>>(
                build_uniform_break_points(theta2_min, theta2_max, theta2_ncells));

        ddc::init_discrete_space<GridR<2>>(SplineInterpPointsR<2>::get_sampling<GridR<2>>());
        ddc::init_discrete_space<GridTheta<2>>(
                SplineInterpPointsTheta<2>::get_sampling<GridTheta<2>>());
    }


    // Initialisation methods --------------------------------------------------------------------
    template <class ExecSpace>
    DConstField<
            typename Patch1::IdxRangeBS12,
            std::experimental::layout_right,
            typename ExecSpace::memory_space>
    set_spline_1(
            DField<typename Patch1::IdxRangeBS12,
                   std::experimental::layout_right,
                   typename ExecSpace::memory_space> const& function_1_coef)
    {
        SplineRThetaBuilder<1, ExecSpace> const builder_1(idx_range_rtheta1);

        // Exact function
        DFieldMem<
                Patch1::IdxRange12,
                ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
                function_1_alloc(idx_range_rtheta1);
        DField<Patch1::IdxRange12,
               std::experimental::layout_right,
               typename ExecSpace::memory_space>
                function_1 = get_field(function_1_alloc);

        ddc::parallel_for_each(
                ExecSpace(),
                idx_range_rtheta1,
                KOKKOS_LAMBDA(Patch1::Idx12 const idx) {
                    double const r = ddc::coordinate(Patch1::Idx1(idx));
                    double const theta = ddc::coordinate(Patch1::Idx2(idx));
                    function_1(idx) = r * Kokkos::sin(theta);
                });

        // Build the spline representations on each patch
        builder_1(function_1_coef, get_const_field(function_1));
        return get_const_field(function_1_coef);
    }

    template <class ExecSpace>
    DConstField<
            typename Patch2::IdxRangeBS12,
            std::experimental::layout_right,
            typename ExecSpace::memory_space>
    set_spline_2(
            DField<typename Patch2::IdxRangeBS12,
                   std::experimental::layout_right,
                   typename ExecSpace::memory_space> const& function_2_coef)
    {
        SplineRThetaBuilder<2, ExecSpace> const builder_2(idx_range_rtheta2);

        // Exact function
        DFieldMem<
                Patch2::IdxRange12,
                ddc::KokkosAllocator<double, typename ExecSpace::memory_space>>
                function_2_alloc(idx_range_rtheta2);
        DField<Patch2::IdxRange12,
               std::experimental::layout_right,
               typename ExecSpace::memory_space>
                function_2 = get_field(function_2_alloc);

        ddc::parallel_for_each(
                ExecSpace(),
                idx_range_rtheta2,
                KOKKOS_LAMBDA(Patch2::Idx12 const idx) {
                    double const r = ddc::coordinate(Patch2::Idx1(idx));
                    double const theta = ddc::coordinate(Patch2::Idx2(idx));
                    function_2(idx) = r * r * Kokkos::sin(theta);
                });

        // Build the spline representations on each patch
        builder_2(function_2_coef, get_const_field(function_2));
        return get_const_field(function_2_coef);
    }
};

template <class PatchLocator>
void test_on_device(
        typename Patch1::Coord1 const& r1_min,
        typename Patch2::Coord1 const& r2_max,
        MultipatchType<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines)
{
    ConstantExtrapolationRuleOnion<PatchLocator> const_extrapol_rule(r1_min, r2_max);

    int constexpr outside_max = PatchLocator::outside_rmax_domain;
    int constexpr outside_min = PatchLocator::outside_rmin_domain;

    ddc::ConstantExtrapolationRule<typename Patch1::Dim1, typename Patch1::Dim2>
            local_extrapol_rule_min(r1_min);
    ddc::ConstantExtrapolationRule<typename Patch2::Dim1, typename Patch2::Dim2>
            local_extrapol_rule_max(r2_max);

    // Outside outter ring.
    typename Patch1::Coord12 outside_coord_max_1(2.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 outside_coord_max_2(2.5, 3 / 2. * M_PI);

    // Inside inner ring
    typename Patch1::Coord12 outside_coord_min_1(.1, 1.6738 * M_PI);
    typename Patch2::Coord12 outside_coord_min_2(.1, 1.6738 * M_PI);

    std::size_t failed_attempt = 0;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<DeviceExecSpace>(0, 1),
            KOKKOS_LAMBDA(int const i, std::size_t& attempt) {
                double expected_value = local_extrapol_rule_max(
                        outside_coord_max_1,
                        splines.template get<Patch2>());
                double extrapol_value
                        = const_extrapol_rule(outside_coord_max_1, splines, outside_max);
                if (expected_value != extrapol_value) {
                    attempt++;
                }

                extrapol_value = const_extrapol_rule(outside_coord_max_2, splines, outside_max);
                if (expected_value != extrapol_value) {
                    attempt++;
                }

                expected_value = local_extrapol_rule_min(
                        outside_coord_min_1,
                        splines.template get<Patch1>());
                extrapol_value = const_extrapol_rule(outside_coord_min_1, splines, outside_min);
                if (expected_value != extrapol_value) {
                    attempt++;
                }

                extrapol_value = const_extrapol_rule(outside_coord_min_2, splines, outside_min);
                if (expected_value != extrapol_value) {
                    attempt++;
                }
            },
            failed_attempt);
    EXPECT_EQ(failed_attempt, 0);
};
} // namespace



TEST_F(MultipatchConstantExtrapolationRuleOnionTest, HostTest)
{
    auto function_1_coef_host = ddc::create_mirror_and_copy(splines.template get<Patch1>());
    auto function_2_coef_host = ddc::create_mirror_and_copy(splines.template get<Patch2>());
    MultipatchType<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const
            splines_host(get_field(function_1_coef_host), get_field(function_2_coef_host));

    ConstantExtrapolationRuleOnion<PatchLocator<HostExecSpace>> const_extrapol_rule(r1_min, r2_max);

    int constexpr outside_max = PatchLocator<HostExecSpace>::outside_rmax_domain;
    int constexpr outside_min = PatchLocator<HostExecSpace>::outside_rmin_domain;

    ddc::ConstantExtrapolationRule<typename Patch1::Dim1, typename Patch1::Dim2>
            local_extrapol_rule_min(r1_min);
    ddc::ConstantExtrapolationRule<typename Patch2::Dim1, typename Patch2::Dim2>
            local_extrapol_rule_max(r2_max);

    double expected_value;
    double extrapol_value;

    // Outside outter ring.
    typename Patch1::Coord12 outside_coord_max_1(2.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 outside_coord_max_2(2.5, 3 / 2. * M_PI);

    expected_value
            = local_extrapol_rule_max(outside_coord_max_1, splines_host.template get<Patch2>());
    extrapol_value = const_extrapol_rule(outside_coord_max_1, splines_host, outside_max);
    EXPECT_EQ(extrapol_value, expected_value);

    extrapol_value = const_extrapol_rule(outside_coord_max_2, splines_host, outside_max);
    EXPECT_EQ(extrapol_value, expected_value);

    // Inside inner ring
    typename Patch1::Coord12 outside_coord_min_1(.1, 1.6738 * M_PI);
    typename Patch2::Coord12 outside_coord_min_2(.1, 1.6738 * M_PI);

    expected_value
            = local_extrapol_rule_min(outside_coord_min_1, splines_host.template get<Patch1>());
    extrapol_value = const_extrapol_rule(outside_coord_min_1, splines_host, outside_min);
    EXPECT_EQ(extrapol_value, expected_value);

    extrapol_value = const_extrapol_rule(outside_coord_min_2, splines_host, outside_min);
    EXPECT_EQ(extrapol_value, expected_value);
}


TEST_F(MultipatchConstantExtrapolationRuleOnionTest, DeviceTest)
{
    test_on_device<PatchLocator<DeviceExecSpace>>(r1_min, r2_max, splines);
}