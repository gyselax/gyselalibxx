// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "kinetic_source.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"
#include "trapezoid_quadrature.hpp"

TEST(KineticSource, Moments)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(100);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IdxStepVx const vx_size(30);

    IdxStepSp const nb_species(2);
    IdxRangeSp const idx_range_sp(IdxSp(0), nb_species);
    IdxSp const my_iion = idx_range_sp.front();
    IdxSp const my_ielec = idx_range_sp.back();

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx gridvx(SplineInterpPointsVx::get_domain<GridVx>());

    SplineXBuilder const builder_x(gridx);
    SplineVxBuilder const builder_vx(gridvx);

    IdxRangeSpXVx const mesh(IdxRangeSp(my_iion, IdxStepSp(1)), gridx, gridvx);

    host_t<DFieldMemX> quadrature_coeffs_x
            = trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(gridx);
    host_t<DFieldMemVx> quadrature_coeffs_vx
            = trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(gridvx);
    host_t<Quadrature<IdxRangeX>> const integrate_x(get_const_field(quadrature_coeffs_x));
    host_t<Quadrature<IdxRangeVx>> const integrate_v(get_const_field(quadrature_coeffs_vx));

    host_t<DFieldMemSp> charges(idx_range_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldMemSp> masses(idx_range_sp);
    ddc::parallel_fill(masses, 1.);

    // Initialisation of the distribution function
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
    DFieldMemSpXVx allfdistribu(mesh);

    // Initialisation of the distribution function
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

    kinetic_source(get_field(allfdistribu), deltat);
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));

    host_t<DFieldMemX> density(gridx);
    host_t<DFieldMemX> fluid_velocity(gridx);
    host_t<DFieldMemX> temperature(gridx);

    host_t<DFieldMemVx> values_density(gridvx);
    host_t<DFieldMemVx> values_fluid_velocity(gridvx);
    host_t<DFieldMemVx> values_temperature(gridvx);
    ddc::for_each(gridx, [&](IdxX const ix) {
        // density
        ddc::parallel_deepcopy(values_density, allfdistribu_host[idx_range_sp.front()][ix]);
        density(ix)
                = integrate_v(Kokkos::DefaultHostExecutionSpace(), get_const_field(values_density));

        // fluid velocity
        ddc::for_each(gridvx, [&](IdxVx const iv) {
            values_fluid_velocity(iv) = values_density(iv) * ddc::coordinate(iv);
        });
        fluid_velocity(ix) = integrate_v(
                                     Kokkos::DefaultHostExecutionSpace(),
                                     get_const_field(values_fluid_velocity))
                             / density(ix);

        // temperature
        ddc::for_each(gridvx, [&](IdxVx const iv) {
            values_temperature(iv)
                    = values_density(iv) * std::pow(ddc::coordinate(iv) - fluid_velocity(ix), 2);
        });
        temperature(ix) = integrate_v(
                                  Kokkos::DefaultHostExecutionSpace(),
                                  get_const_field(values_temperature))
                          / density(ix);
    });

    // source amplitude
    double error_source_amplitude
            = integrate_x(Kokkos::DefaultHostExecutionSpace(), get_const_field(density))
              - source_amplitude;

    double error_fluid_velocity(0);
    double error_temperature(0);
    ddc::for_each(gridx, [&](IdxX const ix) {
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
