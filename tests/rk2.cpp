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
#include "quadrature.hpp"
#include "rk2_solver.hpp"
#include "trapezoid_quadrature.hpp"

TEST(Rk2, Rk2_uniform_vx)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_size(100);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IVectVx const vx_size(30);

    IVectSp const nb_kinspecies(1);

    IDomainSp const dom_sp(IndexSp(0), nb_kinspecies);

    // Creating mesh & supports
    init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(InterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(InterpPointsVx::get_sampling());

    IDomainX gridx(InterpPointsX::get_domain());
    IDomainVx gridvx(InterpPointsVx::get_domain());

    SplineXBuilder const builder_x(gridx);
    SplineVxBuilder const builder_vx(gridvx);

    IDomainSp const gridsp = dom_sp;

    Quadrature<IDimX> const integrate_x(trapezoid_quadrature_coefficients(gridx));
    Quadrature<IDimVx> const integrate_v(trapezoid_quadrature_coefficients(gridvx));

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    SplineEvaluator<BSplinesX> const
            spline_x_evaluator(g_null_boundary<BSplinesX>, g_null_boundary<BSplinesX>);

    SplineEvaluator<BSplinesVx> const
            spline_vx_evaluator(g_null_boundary<BSplinesVx>, g_null_boundary<BSplinesVx>);

    // Initialization of the distribution function --> fill values
    DFieldSpXVx allfdistribu(mesh);
    ddc::fill(allfdistribu, 0.);

    double const tolerance = 1e-20;
    double const deltat = 0.1;

    DFieldVx df(gridvx);
    RK2_solver solver = RK2_solver(
            [](DSpanVx df, DViewSpXVx allfdistribu, double const time, IndexSpX const ispx) {
                ddc::fill(df, time);
            });
    solver(allfdistribu, deltat);

    double const prediction = 0.5 * deltat * deltat;

    for_each(policies::parallel_host, allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        double const val = std::fabs(allfdistribu(ispxvx) - prediction);
        EXPECT_LE(val, tolerance);
    });
}
