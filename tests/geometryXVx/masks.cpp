// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.hpp"
#include "mask_tanh.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

TEST(Masks, Ordering)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());

    IDomainX gridx(SplineInterpPointsX::get_domain());

    SplineXBuilder_1d const builder_x(gridx);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<RDimX> extrapolation_rule_min;
    ddc::PeriodicExtrapolationRule<RDimX> extrapolation_rule_max;
#else
    ddc::ConstantExtrapolationRule<RDimX> extrapolation_rule_min(x_min);
    ddc::ConstantExtrapolationRule<RDimX> extrapolation_rule_max(x_max);
#endif
    SplineXEvaluator_1d const spline_x_evaluator(extrapolation_rule_min, extrapolation_rule_max);

    double const extent = 0.25;
    double const stiffness = 1e-2;

    IndexX const middle(gridx.size() / 2);
    IndexX const tenth(gridx.size() / 10);

    double const tolerance = 1e-10;
    // tests if mask is one inside [x_left, x_right], zero outside
    // with x_left = x_min + Lx*extent
    //      x_right = x_min + Lx - Lx*extent
    //      Lx = gridx total length
    host_t<DFieldX> mask = mask_tanh(gridx, extent, stiffness, MaskType::Normal, false);
    EXPECT_LE(std::fabs(mask(middle) - 1.0), tolerance);
    EXPECT_LE(std::fabs(mask(tenth)), tolerance);

    // same but with a mask that is zero inside [x_left, x_right]
    host_t<DFieldX> mask_inverted = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, false);
    EXPECT_LE(std::fabs(mask_inverted(middle)), tolerance);
    EXPECT_LE(std::fabs(mask_inverted(tenth) - 1.0), tolerance);

    // tests if integral of normalized mask equals 1
    host_t<DFieldX> const quadrature_coeffs = trapezoid_quadrature_coefficients(gridx);
    Quadrature<IDimX> const integrate_x(quadrature_coeffs);

    host_t<DFieldX> mask_normalized = mask_tanh(gridx, extent, stiffness, MaskType::Normal, true);
    EXPECT_LE(std::fabs(integrate_x(mask_normalized) - 1.0), tolerance);
    host_t<DFieldX> mask_normalized_inverted
            = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, true);
    EXPECT_LE(std::fabs(integrate_x(mask_normalized_inverted) - 1.0), tolerance);
}
