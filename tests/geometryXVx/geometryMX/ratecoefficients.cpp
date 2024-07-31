// SPDX-License-Identifier: MIT

#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>

#include "charge_exchange.hpp"
#include "constantfluidinitialization.hpp"
#include "diffusiveneutralsolver.hpp"
#include "geometry.hpp"
#include "ionization.hpp"
#include "maxwellianequilibrium.hpp"
#include "quadrature.hpp"
#include "recombination.hpp"
#include "species_info.hpp"


/**
 * This test initializes a neutral density with a flat spatial profile 
 * and constant quantities for the kinetic plasma species densities, temperature, etc.
 * Then the time derivative of the diffusive neutral model is computed using the solver, 
 * and analytically. The two expressions for the derivative are compared and the relative 
 * difference between the two is checked to be within a percent.
*/
static void TestDiffusiveNeutralsRateCoefficients()
{
    CoordX const x_min(0.0);
    CoordX const x_max(10);
    IdxStepX const x_size(512);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IdxStepVx const vx_size(50);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX meshX(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx meshVx(SplineInterpPointsVx::get_domain<GridVx>());
    IdxRangeXVx meshXVx(meshX, meshVx);

    // Kinetic and neutral species index range initialization
    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const dom_kinsp(IdxSp(0), nb_kinspecies);

    IdxStepSp const nb_fluidspecies(1);
    IdxRangeSp const dom_fluidsp(IdxSp(dom_kinsp.back() + 1), nb_fluidspecies);

    IdxRangeSp const dom_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    host_t<DFieldMemSp> masses(dom_allsp);
    host_t<DFieldSp> kinetic_masses = masses[dom_kinsp];
    host_t<DFieldSp> fluid_masses = masses[dom_fluidsp];

    host_t<DFieldMemSp> charges(dom_allsp);
    host_t<DFieldSp> kinetic_charges = charges[dom_kinsp];
    host_t<DFieldSp> fluid_charges = charges[dom_fluidsp];

    IdxSp const my_iion = dom_kinsp.front();
    IdxSp const my_ielec = dom_kinsp.back();
    IdxSp const my_ifluid = dom_fluidsp.front();

    kinetic_charges(my_ielec) = -1.;
    kinetic_charges(my_iion) = 1.;

    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(my_ielec) = mass_elec;
    kinetic_masses(my_iion) = mass_ion;

    // neutrals charge is zero
    fluid_charges(my_ifluid) = 0.;

    double const neutral_mass(1.);
    fluid_masses(my_ifluid) = neutral_mass;

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    // Moments index range initialization
    IdxStepMom const nb_fluid_moments(1);
    IdxRangeMom const meshM(IdxMom(0), nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();

    IdxRangeSpX dom_fluidspx = IdxRangeSpX(dom_fluidsp, meshX);

    ChargeExchangeRate charge_exchange(1.);
    IonizationRate ionization(1.);
    RecombinationRate recombination(1.);

    DFieldMemSpMomX neutrals_alloc(IdxRangeSpMomX(dom_fluidsp, meshM, meshX));
    DFieldSpMomX neutrals = get_field(neutrals_alloc);

    host_t<DFieldMemSpMom> moments_init(IdxRangeSpMom(dom_fluidsp, meshM));
    ddc::parallel_fill(moments_init, 1.);
    ConstantFluidInitialization fluid_init(moments_init);
    fluid_init(neutrals);

    DFieldMemSpMomX derivative_alloc(get_idx_range(neutrals));
    DFieldSpMomX derivative = get_field(derivative_alloc);

    DFieldMemSpX kinsp_density_alloc(IdxRangeSpX(dom_kinsp, meshX));
    DFieldMemSpX kinsp_velocity_alloc(IdxRangeSpX(dom_kinsp, meshX));
    DFieldMemSpX kinsp_temperature_alloc(IdxRangeSpX(dom_kinsp, meshX));

    DFieldSpX kinsp_density = get_field(kinsp_density_alloc);
    DFieldSpX kinsp_velocity = get_field(kinsp_velocity_alloc);
    DFieldSpX kinsp_temperature = get_field(kinsp_temperature_alloc);

    double const kinsp_density_eq(1.);
    double const kinsp_velocity_eq(0.0);
    double const kinsp_temperature_eq(1.);
    ddc::parallel_fill(kinsp_density, kinsp_density_eq);
    ddc::parallel_fill(kinsp_velocity, kinsp_velocity_eq);
    ddc::parallel_fill(kinsp_temperature, kinsp_temperature_eq);

    // building reaction rates
    DFieldMemSpX charge_exchange_rate_alloc(dom_fluidspx);
    DFieldMemSpX ionization_rate_alloc(dom_fluidspx);
    DFieldMemSpX recombination_rate_alloc(dom_fluidspx);

    DFieldSpX charge_exchange_rate = get_field(charge_exchange_rate_alloc);
    DFieldSpX ionization_rate = get_field(ionization_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);

    charge_exchange(charge_exchange_rate, kinsp_density, kinsp_temperature);
    ionization(ionization_rate, kinsp_density, kinsp_temperature);
    recombination(recombination_rate, kinsp_density, kinsp_temperature);

    double mean_cx_rate = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            0.,
            ddc::reducer::sum<double>(),
            charge_exchange_rate);

    double mean_i_rate = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            0.,
            ddc::reducer::sum<double>(),
            ionization_rate);

    double mean_r_rate = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            0.,
            ddc::reducer::sum<double>(),
            recombination_rate);

    mean_cx_rate /= meshX.size();
    mean_i_rate /= meshX.size();
    mean_r_rate /= meshX.size();

    EXPECT_NEAR(mean_cx_rate, 1.783406341061044, 1e-13);
    EXPECT_NEAR(mean_i_rate, 1.130359390036803, 1e-13);
    EXPECT_NEAR(mean_r_rate, 7.638123065868132e-06, 1e-13);

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}

TEST(GeometryMX, NeutralsRateCoefficients)
{
    TestDiffusiveNeutralsRateCoefficients();
}
