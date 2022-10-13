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
#include "irighthandside.hpp"
#include "kinetic_source.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"
#include "splitrighthandsidesolver.hpp"
#include "trapezoid_quadrature.hpp"

TEST(KineticSource, Ordering)
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

    SplineXBuilder const builder_x;

    SplineVxBuilder const builder_vx;

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSp const gridsp = dom_sp;

    Quadrature<IDimX> const integrate_x(trapezoid_quadrature_coefficients(gridx));
    Quadrature<IDimVx> const integrate_v(trapezoid_quadrature_coefficients(gridvx));

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    SplineEvaluator<BSplinesX> const
            spline_x_evaluator(g_null_boundary<BSplinesX>, g_null_boundary<BSplinesX>);

    SplineEvaluator<BSplinesVx> const
            spline_vx_evaluator(g_null_boundary<BSplinesVx>, g_null_boundary<BSplinesVx>);

    FieldSp<int> charges(dom_sp);
    charges(dom_sp.front()) = 1;
    DFieldSp masses(dom_sp);
    masses(dom_sp.front()) = 1.0;
    FieldSp<int> init_perturb_mode(dom_sp);
    init_perturb_mode(dom_sp.front()) = 0;
    DFieldSp init_perturb_amplitude(dom_sp);
    init_perturb_amplitude(dom_sp.front()) = 0.0;

    // Initialization of the distribution function
    init_discrete_space<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));

    DFieldSpXVx allfdistribu(mesh);

    // Initialization of the distribution function --> fill values
    double fdistribu_val = 0.;
    for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        allfdistribu(ispxvx) = fdistribu_val;
    });

    // Maxwellian source test
    double const px_source = 0.2;
    double const dx_source = 0.1;
    double const source_amplitude = 1.;
    double const density_amplitude = 1;
    double const energy_amplitude = 1;
    double const temperature_source = 0.5;
    //
    // --> Algorithm info
    double const deltat = 1.;

    Kinetic_source const kinetic_source(
            gridx,
            gridvx,
            px_source,
            dx_source,
            source_amplitude,
            density_amplitude,
            energy_amplitude,
            temperature_source);

    kinetic_source(allfdistribu, deltat);

    DFieldX density(gridx);
    DFieldX fluid_velocity(gridx);
    DFieldX temperature(gridx);

    DFieldVx values_density(gridvx);
    DFieldVx values_fluid_velocity(gridvx);
    DFieldVx values_temperature(gridvx);
    for_each(gridx, [&](IndexX const ix) {
        // density
        deepcopy(values_density, allfdistribu[dom_sp.front()][ix]);
        density(ix) = integrate_v(values_density);

        // fluid velocity
        for_each(gridvx, [&](IndexVx const iv) {
            values_fluid_velocity(iv) = values_density(iv) * coordinate(iv);
        });
        fluid_velocity(ix) = integrate_v(values_fluid_velocity) / density(ix);

        // temperature
        for_each(gridvx, [&](IndexVx const iv) {
            values_temperature(iv)
                    = values_density(iv) * std::pow(coordinate(iv) - fluid_velocity(ix), 2);
        });
        temperature(ix) = integrate_v(values_temperature) / density(ix);
    });

    // source amplitude
    double error_source_amplitude = integrate_x(density) - source_amplitude;

    double error_fluid_velocity(0);
    double error_temperature(0);
    for_each(gridx, [&](IndexX const ix) {
        error_fluid_velocity = std::fmax(std::fabs(fluid_velocity(ix)), error_fluid_velocity);
        error_temperature
                = std::fmax(std::fabs(temperature(ix) - temperature_source), error_temperature);
    });
    EXPECT_LE(error_source_amplitude, 1e-3);
    EXPECT_LE(error_fluid_velocity, 1e-8);
    EXPECT_LE(error_temperature, 1e-8);
}
