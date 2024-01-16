// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <paraconf.h>
#include <pdi.h>

#include "chargedensitycalculator.hpp"
#include "fftpoissonsolver.hpp"
#include "geometry.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"

TEST(FftPoissonSolver, CosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IVectX const x_size(100);

    CoordVx const vx_min(-0.5);
    CoordVx const vx_max(0.5);
    IVectVx const vx_size(11);

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

    SplineBuilder<BSplinesX, IDimX, BoundCond::PERIODIC, BoundCond::PERIODIC> const builder_x(
            interpolation_domain_x);

    SplineVxBuilder const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSp const gridsp = IDomainSp(my_iion, IVectSp(1));

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    ddc::init_fourier_space<RDimX>(ddc::select<IDimX>(mesh));

    // Creating operators

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

    DFieldVx const quadrature_coeffs = neumann_spline_quadrature_coefficients(gridvx, builder_vx);
    Quadrature<IDimVx> const integrate_v(quadrature_coeffs);
    ChargeDensityCalculator rhs(integrate_v);
    FftPoissonSolver poisson(builder_x, spline_x_evaluator, rhs);

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
    device_t<DFieldX> electrostatic_potential_device(gridx);
    device_t<DFieldX> electric_field_device(gridx);
    device_t<DFieldSpXVx> allfdistribu_device(mesh);

    ddc::deepcopy(allfdistribu_device, allfdistribu);
    poisson(electrostatic_potential_device, electric_field_device, allfdistribu_device);
    ddc::deepcopy(electric_field, electric_field_device);
    ddc::deepcopy(electrostatic_potential, electrostatic_potential_device);

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