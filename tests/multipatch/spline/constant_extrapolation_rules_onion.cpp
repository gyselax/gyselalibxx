// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_non_uniform.hpp"
#include "constant_extrapolation_rules_onion.hpp"
#include "ddc_helper.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "onion_patch_locator.hpp"
#include "physical_geometry.hpp"
#include "spline_testing_tools.hpp"
#include "types.hpp"


namespace {
template <class PatchLocator>
void test_on_device(
        typename Patch1::Coord1 const& r1_min,
        typename Patch2::Coord1 const& r2_max,
        MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines)
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



TEST_F(MultipatchSplineOnionShapeTest, ConstantExtrapolationRuleOnionShapeHostTest)
{
    auto function_1_coef_host = ddc::create_mirror_and_copy(splines.template get<Patch1>());
    auto function_2_coef_host = ddc::create_mirror_and_copy(splines.template get<Patch2>());
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const
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


TEST_F(MultipatchSplineOnionShapeTest, ConstantExtrapolationRuleOnionShapeDeviceTest)
{
    test_on_device<PatchLocator<DeviceExecSpace>>(r1_min, r2_max, splines);
}
