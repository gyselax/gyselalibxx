// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <paraconf.h>
#include <pdi.h>

#include "femperiodicpoissonsolver.hpp"
#include "geometry.hpp"
#include "species_info.hpp"

TEST(FemPeriodicPoissonSolver, CosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IVectX const x_size(100);

    CoordVx const vx_min(-0.5);
    CoordVx const vx_max(0.5);
    IVectVx const vx_size(10);

    IVectSp const nb_kinspecies(1);

    IDomainSp const dom_sp(nb_kinspecies);

    // Creating mesh & supports
    init_discretization<BSplinesX>(x_min, x_max, x_size);

    init_discretization<BSplinesVx>(vx_min, vx_max, vx_size);

    SplineXBuilder const builder_x;

    SplineVxBuilder const builder_vx;

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSp const gridsp = dom_sp;

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    SplineEvaluator<BSplinesX> const spline_x_evaluator(
            NullBoundaryValue<BSplinesX>::value,
            NullBoundaryValue<BSplinesX>::value);

    SplineEvaluator<BSplinesVx> const spline_vx_evaluator(
            NullBoundaryValue<BSplinesVx>::value,
            NullBoundaryValue<BSplinesVx>::value);

    FieldSp<int> charges(dom_sp);
    charges(dom_sp.front()) = 1;
    DFieldSp masses(dom_sp);
    masses(dom_sp.front()) = 1.0;
    FieldSp<int> init_perturb_mode(dom_sp);
    init_perturb_mode(dom_sp.front()) = 0;
    DFieldSp init_perturb_amplitude(dom_sp);
    init_perturb_amplitude(dom_sp.front()) = 0.0;

    // Initialization of the distribution function
    init_discretization<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));

    FemPeriodicPoissonSolver
            poisson(builder_x, spline_x_evaluator, builder_vx, spline_vx_evaluator);

    DFieldX electrostatic_potential(gridx);
    DFieldX electric_field(gridx);
    DFieldSpXVx allfdistribu(mesh);

    // Initialization of the distribution function --> fill values
    for (IndexSp const isp : gridsp) {
        for (IndexX const ix : gridx) {
            double fdistribu_val = cos(to_real(ix));
            for (IndexVx const iv : gridvx) {
                allfdistribu(isp, ix, iv) = fdistribu_val;
            }
        }
    }

    poisson(electrostatic_potential, electric_field, allfdistribu);

    double error_pot = 0.0;
    double error_field = 0.0;

    for (IndexX const ix : gridx) {
        double const exact_pot = cos(to_real(ix));
        error_pot = fmax(fabs(electrostatic_potential(ix) - exact_pot), error_pot);
        double const exact_field = sin(to_real(ix));
        error_field = fmax(fabs(electric_field(ix) - exact_field), error_field);
    }
    EXPECT_LE(error_pot, 1e-8);
    EXPECT_LE(error_field, 1e-6);
}
