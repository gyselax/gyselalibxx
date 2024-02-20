// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <paraconf.h>
#include <pdi.h>

#include "chargedensitycalculator.hpp"
#include "femnonperiodicpoissonsolver.hpp"
#include "geometry.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"

TEST(FemNonPeriodicPoissonSolver, Ordering)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
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

    SplineXBuilder_1d const builder_x(interpolation_domain_x);

    SplineVxBuilder_1d const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSp const gridsp = IDomainSp(my_iion, IVectSp(1));

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    SplineXEvaluator_1d const spline_x_evaluator(
            builder_x.spline_domain(),
            ddc::ConstantExtrapolationRule<RDimX>(x_min),
            ddc::ConstantExtrapolationRule<RDimX>(x_max));

    SplineVxEvaluator_1d const spline_vx_evaluator(
            builder_vx.spline_domain(),
            ddc::ConstantExtrapolationRule<RDimVx>(vx_min),
            ddc::ConstantExtrapolationRule<RDimVx>(vx_min));

    FieldSp<int> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    DFieldSp masses(dom_sp);
    ddc::fill(masses, 1);

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

    DFieldVx const quadrature_coeffs = neumann_spline_quadrature_coefficients(gridvx, builder_vx);
    Quadrature<IDimVx> const integrate_v(quadrature_coeffs);
    ChargeDensityCalculator rhs(integrate_v);
    FemNonPeriodicPoissonSolver poisson(builder_x, spline_x_evaluator, rhs);

    DFieldX electrostatic_potential(gridx);
    DFieldX electric_field(gridx);
    DFieldSpXVx allfdistribu(mesh);

    // Initialization of the distribution function --> fill values
    for (IndexSp const isp : gridsp) {
        for (IndexX const ix : gridx) {
            double fdistribu_val = sin(ddc::coordinate(ix));
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
        double const exact_pot = sin(ddc::coordinate(ix));
        error_pot = fmax(fabs(electrostatic_potential(ix) - exact_pot), error_pot);
        double const exact_field = -cos(ddc::coordinate(ix));
        error_field = fmax(fabs(electric_field(ix) - exact_field), error_field);
    }
    EXPECT_LE(error_pot, 1e-2);
    EXPECT_LE(error_field, 1e-1);
}
