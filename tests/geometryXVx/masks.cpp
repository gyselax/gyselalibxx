// SPDX-License-Identifier: MIT
#include <cmath>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "mask_tanh.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

TEST(Masks, Ordering)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());

    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());

    SplineXBuilder const builder_x(gridx);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<X> extrapolation_rule_min;
    ddc::PeriodicExtrapolationRule<X> extrapolation_rule_max;
#else
    ddc::ConstantExtrapolationRule<X> extrapolation_rule_min(x_min);
    ddc::ConstantExtrapolationRule<X> extrapolation_rule_max(x_max);
#endif
    SplineXEvaluator const spline_x_evaluator(extrapolation_rule_min, extrapolation_rule_max);

    double const extent = 0.25;
    double const stiffness = 1e-2;

    IdxX const middle(gridx.size() / 2);
    IdxX const tenth(gridx.size() / 10);

    double const tolerance = 1e-10;
    // tests if mask is one inside [x_left, x_right], zero outside
    // with x_left = x_min + Lx*extent
    //      x_right = x_min + Lx - Lx*extent
    //      Lx = gridx total length
    host_t<DFieldMemX> mask = mask_tanh(gridx, extent, stiffness, MaskType::Normal, false);
    EXPECT_LE(std::fabs(mask(middle) - 1.0), tolerance);
    EXPECT_LE(std::fabs(mask(tenth)), tolerance);

    // same but with a mask that is zero inside [x_left, x_right]
    host_t<DFieldMemX> mask_inverted
            = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, false);
    EXPECT_LE(std::fabs(mask_inverted(middle)), tolerance);
    EXPECT_LE(std::fabs(mask_inverted(tenth) - 1.0), tolerance);

    // tests if integral of normalised mask equals 1
    DFieldMemX const quadrature_coeffs(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridx));
    Quadrature<IdxRangeX> const integrate_x(get_const_field(quadrature_coeffs));

    host_t<DFieldMemX> mask_normalised_host
            = mask_tanh(gridx, extent, stiffness, MaskType::Normal, true);
    auto mask_normalised = ddc::create_mirror_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(mask_normalised_host));
    double const mask_integrated
            = integrate_x(Kokkos::DefaultExecutionSpace(), get_const_field(mask_normalised));
    EXPECT_LE(std::fabs(mask_integrated - 1.0), tolerance);
    host_t<DFieldMemX> mask_normalised_inverted_host
            = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, true);
    auto mask_normalised_inverted = ddc::create_mirror_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(mask_normalised_inverted_host));
    double const mask_inverted_integrated = integrate_x(
            Kokkos::DefaultExecutionSpace(),
            get_const_field(mask_normalised_inverted));
    EXPECT_LE(std::fabs(mask_inverted_integrated - 1.0), tolerance);
}
