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

    IVectSp const nb_species(2);
    IDomainSp const dom_sp(IndexSp(0), nb_species);
    IndexSp const my_ielec = dom_sp.front();
    IndexSp const my_iion = dom_sp.back();

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());
    ddc::DiscreteDomain<IDimX> interpolation_domain_x(SplineInterpPointsX::get_domain());
    ddc::DiscreteDomain<IDimVx> interpolation_domain_vx(SplineInterpPointsVx::get_domain());

    SplineXBuilder const builder_x(interpolation_domain_x);

    SplineVxBuilder const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSp const gridsp = IDomainSp(my_iion, IVectSp(1));

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    SplineEvaluator<BSplinesX> const
            spline_x_evaluator(g_null_boundary<BSplinesX>, g_null_boundary<BSplinesX>);

    SplineEvaluator<BSplinesVx> const
            spline_vx_evaluator(g_null_boundary<BSplinesVx>, g_null_boundary<BSplinesVx>);

    FieldSp<int> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    DFieldSp masses(dom_sp);
    ddc::fill(masses, 1);
    FieldSp<int> init_perturb_mode(dom_sp);
    ddc::fill(init_perturb_mode, 0);
    DFieldSp init_perturb_amplitude(dom_sp);
    ddc::fill(init_perturb_amplitude, 0);

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(
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
            double fdistribu_val = cos(ddc::coordinate(ix));
            for (IndexVx const iv : gridvx) {
                allfdistribu(isp, ix, iv) = fdistribu_val;
            }
        }
    }

    poisson(electrostatic_potential, electric_field, allfdistribu);

    double error_pot = 0.0;
    double error_field = 0.0;

    for (IndexX const ix : gridx) {
        double const exact_pot = cos(ddc::coordinate(ix));
        error_pot = fmax(fabs(electrostatic_potential(ix) - exact_pot), error_pot);
        double const exact_field = sin(ddc::coordinate(ix));
        error_field = fmax(fabs(electric_field(ix) - exact_field), error_field);
    }
    EXPECT_LE(error_pot, 1e-8);
    EXPECT_LE(error_field, 1e-6);
}
