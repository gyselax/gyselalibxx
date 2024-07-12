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
    IVectX const x_size(512);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IVectVx const vx_size(50);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());

    IDomainX meshX(SplineInterpPointsX::get_domain<IDimX>());
    IDomainVx meshVx(SplineInterpPointsVx::get_domain<IDimVx>());
    IDomainXVx meshXVx(meshX, meshVx);

    // Kinetic and neutral species domain initialization
    IVectSp const nb_kinspecies(2);
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IVectSp const nb_fluidspecies(1);
    IDomainSp const dom_fluidsp(IndexSp(dom_kinsp.back() + 1), nb_fluidspecies);

    IDomainSp const dom_allsp(IndexSp(0), nb_kinspecies + nb_fluidspecies);

    host_t<DFieldSp> masses(dom_allsp);
    host_t<DSpanSp> kinetic_masses = masses[dom_kinsp];
    host_t<DSpanSp> fluid_masses = masses[dom_fluidsp];

    host_t<FieldSp<int>> charges(dom_allsp);
    host_t<SpanSp<int>> kinetic_charges = charges[dom_kinsp];
    host_t<SpanSp<int>> fluid_charges = charges[dom_fluidsp];

    IndexSp const my_iion = dom_kinsp.front();
    IndexSp const my_ielec = dom_kinsp.back();
    IndexSp const my_ifluid = dom_fluidsp.front();

    kinetic_charges(my_ielec) = -1;
    kinetic_charges(my_iion) = 1;

    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(my_ielec) = mass_elec;
    kinetic_masses(my_iion) = mass_ion;

    // neutrals charge is zero
    fluid_charges(my_ifluid) = 0;

    double const neutral_mass(1.);
    fluid_masses(my_ifluid) = neutral_mass;

    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

    // Moments domain initialization
    IVectM const nb_fluid_moments(1);
    IDomainM const meshM(IndexM(0), nb_fluid_moments);
    ddc::init_discrete_space<IDimM>();

    IDomainSpX dom_fluidspx = IDomainSpX(dom_fluidsp, meshX);

    ChargeExchangeRate charge_exchange(1.);
    IonizationRate ionization(1.);
    RecombinationRate recombination(1.);

    DFieldSpMX neutrals_alloc(IDomainSpMX(dom_fluidsp, meshM, meshX));
    DSpanSpMX neutrals = neutrals_alloc.span_view();

    host_t<DFieldSpM> moments_init(IDomainSpM(dom_fluidsp, meshM));
    ddc::parallel_fill(moments_init, 1.);
    ConstantFluidInitialization fluid_init(moments_init);
    fluid_init(neutrals);

    DFieldSpMX derivative_alloc(neutrals.domain());
    DSpanSpMX derivative = derivative_alloc.span_view();

    DFieldSpX kinsp_density_alloc(IDomainSpX(dom_kinsp, meshX));
    DFieldSpX kinsp_velocity_alloc(IDomainSpX(dom_kinsp, meshX));
    DFieldSpX kinsp_temperature_alloc(IDomainSpX(dom_kinsp, meshX));

    DSpanSpX kinsp_density = kinsp_density_alloc.span_view();
    DSpanSpX kinsp_velocity = kinsp_velocity_alloc.span_view();
    DSpanSpX kinsp_temperature = kinsp_temperature_alloc.span_view();

    double const kinsp_density_eq(1.);
    double const kinsp_velocity_eq(0.0);
    double const kinsp_temperature_eq(1.);
    ddc::parallel_fill(kinsp_density, kinsp_density_eq);
    ddc::parallel_fill(kinsp_velocity, kinsp_velocity_eq);
    ddc::parallel_fill(kinsp_temperature, kinsp_temperature_eq);

    // building reaction rates
    DFieldSpX charge_exchange_rate_alloc(dom_fluidspx);
    DFieldSpX ionization_rate_alloc(dom_fluidspx);
    DFieldSpX recombination_rate_alloc(dom_fluidspx);

    DSpanSpX charge_exchange_rate = charge_exchange_rate_alloc.span_view();
    DSpanSpX ionization_rate = ionization_rate_alloc.span_view();
    DSpanSpX recombination_rate = recombination_rate_alloc.span_view();

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
