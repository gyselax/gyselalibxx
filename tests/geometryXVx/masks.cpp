// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

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

    ddc::init_discrete_space<IDimX>(InterpPointsX::get_sampling());

    IDomainX gridx(InterpPointsX::get_domain());

    SplineXBuilder const builder_x(gridx);

    SplineEvaluator<BSplinesX> const
            spline_x_evaluator(g_null_boundary<BSplinesX>, g_null_boundary<BSplinesX>);

    double const extent = 0.25;
    double const stiffness = 1e-2;

    IndexX const middle(gridx.size() / 2);
    IndexX const tenth(gridx.size() / 10);

    double const tolerance = 1e-10;
    // tests if mask is one inside [x_left, x_right], zero outside
    // with x_left = x_min + Lx*extent
    //      x_right = x_min + Lx - Lx*extent
    //      Lx = gridx total length
    DFieldX mask = mask_tanh(gridx, extent, stiffness, MaskType::Normal, false);
    EXPECT_LE(std::fabs(mask(middle) - 1.0), tolerance);
    EXPECT_LE(std::fabs(mask(tenth)), tolerance);

    // same but with a mask that is zero inside [x_left, x_right]
    DFieldX mask_inverted = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, false);
    EXPECT_LE(std::fabs(mask_inverted(middle)), tolerance);
    EXPECT_LE(std::fabs(mask_inverted(tenth) - 1.0), tolerance);

    // tests if integral of normalized mask equals 1
    Quadrature<IDimX> const integrate_x(trapezoid_quadrature_coefficients(gridx));

    DFieldX mask_normalized = mask_tanh(gridx, extent, stiffness, MaskType::Normal, true);
    EXPECT_LE(std::fabs(integrate_x(mask_normalized) - 1.0), tolerance);
    DFieldX mask_normalized_inverted
            = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, true);
    EXPECT_LE(std::fabs(integrate_x(mask_normalized_inverted) - 1.0), tolerance);
}
