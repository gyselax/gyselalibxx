// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <geometry.hpp>
#include <irighthandside.hpp>
#include <kinetic_source.hpp>
#include <pdi.h>
#include <quadrature.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

TEST(KineticSource, Moments)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_size(100);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IVectVx const vx_size(30);

    IVectSp const nb_species(2);
    IDomainSp const dom_sp(IndexSp(0), nb_species);
    IndexSp const my_iion = dom_sp.front();
    IndexSp const my_ielec = dom_sp.back();

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());

    IDomainX interpolation_domain_x(SplineInterpPointsX::get_domain<IDimX>());
    IDomainVx interpolation_domain_vx(SplineInterpPointsVx::get_domain<IDimVx>());

    SplineXBuilder_1d const builder_x(interpolation_domain_x);

    SplineVxBuilder_1d const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSpXVx const mesh(IDomainSp(my_iion, IVectSp(1)), gridx, gridvx);

    host_t<DFieldX> quadrature_coeffs_x = trapezoid_quadrature_coefficients(gridx);
    host_t<DFieldVx> quadrature_coeffs_vx = trapezoid_quadrature_coefficients(gridvx);
    Quadrature<IDimX> const integrate_x(quadrature_coeffs_x);
    Quadrature<IDimVx> const integrate_v(quadrature_coeffs_vx);

    host_t<DFieldSp> charges(dom_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldSp> masses(dom_sp);
    ddc::parallel_fill(masses, 1.);

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));
    DFieldSpXVx allfdistribu(mesh);

    // Initialization of the distribution function
    ddc::parallel_fill(allfdistribu, 0.);

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

    KineticSource const kinetic_source(
            gridx,
            gridvx,
            px_source,
            dx_source,
            source_amplitude,
            density_amplitude,
            energy_amplitude,
            temperature_source);

    kinetic_source(allfdistribu, deltat);
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu.span_view());

    host_t<DFieldX> density(gridx);
    host_t<DFieldX> fluid_velocity(gridx);
    host_t<DFieldX> temperature(gridx);

    host_t<DFieldVx> values_density(gridvx);
    host_t<DFieldVx> values_fluid_velocity(gridvx);
    host_t<DFieldVx> values_temperature(gridvx);
    ddc::for_each(gridx, [&](IndexX const ix) {
        // density
        ddc::parallel_deepcopy(values_density, allfdistribu_host[dom_sp.front()][ix]);
        density(ix) = integrate_v(values_density);

        // fluid velocity
        ddc::for_each(gridvx, [&](IndexVx const iv) {
            values_fluid_velocity(iv) = values_density(iv) * ddc::coordinate(iv);
        });
        fluid_velocity(ix) = integrate_v(values_fluid_velocity) / density(ix);

        // temperature
        ddc::for_each(gridvx, [&](IndexVx const iv) {
            values_temperature(iv)
                    = values_density(iv) * std::pow(ddc::coordinate(iv) - fluid_velocity(ix), 2);
        });
        temperature(ix) = integrate_v(values_temperature) / density(ix);
    });

    // source amplitude
    double error_source_amplitude = integrate_x(density) - source_amplitude;

    double error_fluid_velocity(0);
    double error_temperature(0);
    ddc::for_each(gridx, [&](IndexX const ix) {
        error_fluid_velocity = std::fmax(std::fabs(fluid_velocity(ix)), error_fluid_velocity);
        error_temperature
                = std::fmax(std::fabs(temperature(ix) - temperature_source), error_temperature);
    });
    EXPECT_LE(error_source_amplitude, 1e-3);
    EXPECT_LE(error_fluid_velocity, 1e-8);
    EXPECT_LE(error_temperature, 1e-8);

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
