// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_non_uniform.hpp"
#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "null_extrapolation_rules.hpp"
#include "spline_testing_tools.hpp"
#include "types.hpp"



namespace {
void test_on_device(MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const& splines)
{
    NullExtrapolationRule null_extrapol_rule;

    int constexpr outside = -1;

    // Outside outer ring.
    typename Patch1::Coord12 outside_coord_1(2.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 outside_coord_2(2.5, 3 / 2. * M_PI);

    std::size_t failed_attempt = 0;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<DeviceExecSpace>(0, 1),
            KOKKOS_LAMBDA(int const i, std::size_t& attempt) {
                if ((null_extrapol_rule(outside_coord_1, splines, outside) != 0)
                    || (null_extrapol_rule(outside_coord_2, splines, outside) != 0)) {
                    attempt++;
                }
            },
            failed_attempt);
    EXPECT_EQ(failed_attempt, 0);
};

} // namespace



TEST_F(MultipatchSplineOnionShapeTest, NullExtrapolationRuleHostTest)
{
    auto function_1_coef_host = ddc::create_mirror_and_copy(splines.template get<Patch1>());
    auto function_2_coef_host = ddc::create_mirror_and_copy(splines.template get<Patch2>());
    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patch1, Patch2> const splines_host(
            get_const_field(function_1_coef_host),
            get_const_field(function_2_coef_host));

    NullExtrapolationRule null_extrapol_rule;

    int constexpr outside = -1;

    // Outside outer ring.
    typename Patch1::Coord12 outside_coord_1(2.5, 3 / 2. * M_PI);
    typename Patch2::Coord12 outside_coord_2(2.5, 3 / 2. * M_PI);

    EXPECT_EQ(null_extrapol_rule(outside_coord_1, splines_host, outside), 0);
    EXPECT_EQ(null_extrapol_rule(outside_coord_2, splines_host, outside), 0);
}


TEST_F(MultipatchSplineOnionShapeTest, NullExtrapolationRuleDeviceTest)
{
    test_on_device(splines);
}
