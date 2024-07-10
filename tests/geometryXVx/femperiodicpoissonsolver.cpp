// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <paraconf.h>
#include <pdi.h>

#include "chargedensitycalculator.hpp"
#include "femperiodicqnsolver.hpp"
#include "geometry.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"

TEST(FemPeriodicQNSolver, CosineSource)
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
    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    ddc::DiscreteDomain<IDimX> interpolation_domain_x(SplineInterpPointsX::get_domain<IDimX>());
    ddc::DiscreteDomain<IDimVx> interpolation_domain_vx(SplineInterpPointsVx::get_domain<IDimVx>());

    SplineXBuilder_1d const builder_x(interpolation_domain_x);
    SplineVxBuilder_1d const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSp const gridsp = IDomainSp(my_iion, IVectSp(1));

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    ddc::PeriodicExtrapolationRule<RDimX> x_extrapolation_rule_min;
    ddc::PeriodicExtrapolationRule<RDimX> x_extrapolation_rule_max;
    ddc::ConstantExtrapolationRule<RDimVx> vx_extrapolation_rule_min(vx_min);
    ddc::ConstantExtrapolationRule<RDimVx> vx_extrapolation_rule_max(vx_max);
    SplineXEvaluator_1d const
            spline_x_evaluator(x_extrapolation_rule_min, x_extrapolation_rule_max);

    SplineVxEvaluator_1d const
            spline_vx_evaluator(vx_extrapolation_rule_min, vx_extrapolation_rule_max);
    host_t<FieldSp<int>> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    host_t<DFieldSp> masses(dom_sp);
    ddc::parallel_fill(masses, 1);

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

    DFieldVx const quadrature_coeffs(neumann_spline_quadrature_coefficients(gridvx, builder_vx));
    ChargeDensityCalculator rhs(quadrature_coeffs);
    FemPeriodicQNSolver poisson(builder_x, spline_x_evaluator, rhs);

    host_t<DFieldX> electrostatic_potential_host(gridx);
    host_t<DFieldX> electric_field_host(gridx);
    host_t<DFieldSpXVx> allfdistribu_host(mesh);

    // Initialization of the distribution function --> fill values
    for (IndexSp const isp : gridsp) {
        for (IndexX const ix : gridx) {
            double fdistribu_val = cos(ddc::coordinate(ix));
            for (IndexVx const iv : gridvx) {
                allfdistribu_host(isp, ix, iv) = fdistribu_val;
            }
        }
    }
    DFieldX electrostatic_potential(gridx);
    DFieldX electric_field(gridx);
    DFieldSpXVx allfdistribu(mesh);

    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    poisson(electrostatic_potential, electric_field, allfdistribu);
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);

    double error_pot = 0.0;
    double error_field = 0.0;


    for (IndexX const ix : gridx) {
        double const exact_pot = cos(ddc::coordinate(ix));
        error_pot = fmax(fabs(electrostatic_potential_host(ix) - exact_pot), error_pot);
        double const exact_field = sin(ddc::coordinate(ix));
        error_field = fmax(fabs(electric_field_host(ix) - exact_field), error_field);
    }
    EXPECT_LE(error_pot, 1e-8);
    EXPECT_LE(error_field, 1e-6);
}
