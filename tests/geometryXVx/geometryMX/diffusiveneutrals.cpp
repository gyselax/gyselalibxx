// SPDX-License-Identifier: MIT

#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

#include "constantfluidinitialization.hpp"
#include "constantrate.hpp"
#include "diffusiveneutralsolver.hpp"
#include "geometry.hpp"
#include "maxwellianequilibrium.hpp"
#include "species_info.hpp"

/**
 * This test initializes a neutral density with an exponential spatial variation 
 * and constant quantities for the kinetic species densities, temperature, etc.
 * Then the time derivative of the diffusive neutral model is computed using the solver, 
 * and analytically. The two expressions for the derivative are compared and the relative 
 * difference between the two is checked to be within a percent. 
*/
TEST(GeometryMX, DiffusiveNeutralsDerivative)
{
    CoordX const x_min(0.0);
    CoordX const x_max(10);
    IdxStepX const x_size(512);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IdxStepVx const vx_size(50);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX meshX(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx meshVx(SplineInterpPointsVx::get_domain<GridVx>());
    IdxRangeXVx meshXVx(meshX, meshVx);

    SplineXBuilder const builder_x(meshXVx);
    SplineVxBuilder const builder_vx(meshXVx);

    // Kinetic species index range initialization
    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const dom_kinsp(IdxSp(0), nb_kinspecies);

    IdxSp const my_iion = dom_kinsp.front();
    IdxSp const my_ielec = dom_kinsp.back();

    host_t<DFieldMemSp> kinetic_charges(dom_kinsp);
    kinetic_charges(my_ielec) = -1.;
    kinetic_charges(my_iion) = 1.;

    host_t<DFieldMemSp> kinetic_masses(dom_kinsp);
    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(my_ielec) = mass_elec;
    kinetic_masses(my_iion) = mass_ion;

    // Neutral species index range initialization
    IdxStepSp const nb_fluidspecies(1);
    IdxRangeSp const dom_fluidsp(IdxSp(dom_kinsp.back() + 1), nb_fluidspecies);
    IdxSp const my_ifluid = dom_fluidsp.front();

    // neutrals charge is zero
    host_t<DFieldMemSp> fluid_charges(dom_fluidsp);
    ddc::parallel_fill(fluid_charges, 0.);

    host_t<DFieldMemSp> fluid_masses(dom_fluidsp);
    double const neutral_mass(1.);
    fluid_masses(my_ifluid) = neutral_mass;

    // Create the index range of kinetic species + fluid species
    IdxRangeSp const dom_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains charges of all species
    host_t<DFieldMemSp> charges(dom_allsp);

    // fill the Field with charges of kinetic species
    for (IdxSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }

    // fill the Field with charges of fluid species
    for (IdxSp isp : dom_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }

    // Create a Field that contains masses of kinetic and fluid species
    host_t<DFieldMemSp> masses(dom_allsp);

    // fill the Field with masses of kinetic species
    for (IdxSp isp : dom_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }

    // fill the Field with masses of fluid species
    for (IdxSp isp : dom_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    // Moments index range initialization
    IdxStepMom const nb_fluid_moments(1);
    IdxRangeMom const meshM(IdxMom(0), nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();

    // Neutral species initialization
    double const charge_exchange_val(0.5);
    double const ionization_val(1.);
    double const recombination_val(2.);
    double const normalization_coeff(1.);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
#else
    ddc::ConstantExtrapolationRule<X> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<X> bv_x_max(x_max);
#endif

    ConstantRate charge_exchange(charge_exchange_val);
    ConstantRate ionization(ionization_val);
    ConstantRate recombination(recombination_val);

    SplineXBuilder_1d const spline_x_builder_neutrals(meshX);
    SplineXEvaluator_1d const spline_x_evaluator_neutrals(bv_x_min, bv_x_max);

    DFieldMemVx quadrature_coeffs_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(meshVx));
    DFieldVx const quadrature_coeffs = get_field(quadrature_coeffs_alloc);

    DiffusiveNeutralSolver const neutralsolver(
            charge_exchange,
            ionization,
            recombination,
            normalization_coeff,
            spline_x_builder_neutrals,
            spline_x_evaluator_neutrals,
            quadrature_coeffs);

    host_t<DFieldMemSpMomX> neutrals_init_host(IdxRangeSpMomX(dom_fluidsp, meshM, meshX));
    ddc::for_each(get_idx_range(neutrals_init_host), [&](IdxSpMomX const ispmx) {
        CoordX coordx(ddc::coordinate(ddc::select<GridX>(ispmx)));
        double const lx_2((x_max + x_min) / 2.);
        neutrals_init_host(ispmx) = std::exp(-0.5 * (coordx - lx_2) * (coordx - lx_2));
    });

    DFieldMemSpMomX neutrals_alloc(get_idx_range(neutrals_init_host));
    DFieldSpMomX neutrals = get_field(neutrals_alloc);
    ddc::parallel_deepcopy(neutrals, neutrals_init_host);

    DFieldMemSpMomX derivative_alloc(get_idx_range(neutrals));
    DFieldSpMomX derivative = get_field(derivative_alloc);

    DFieldMemSpX kinsp_density_alloc(IdxRangeSpX(dom_kinsp, meshX));
    DFieldMemSpX kinsp_velocity_alloc(IdxRangeSpX(dom_kinsp, meshX));
    DFieldMemSpX kinsp_temperature_alloc(IdxRangeSpX(dom_kinsp, meshX));

    DFieldSpX kinsp_density = get_field(kinsp_density_alloc);
    DFieldSpX kinsp_velocity = get_field(kinsp_velocity_alloc);
    DFieldSpX kinsp_temperature = get_field(kinsp_temperature_alloc);

    double const kinsp_density_eq(1.);
    double const kinsp_velocity_eq(0.5);
    double const kinsp_temperature_eq(1.);
    ddc::parallel_fill(kinsp_density, kinsp_density_eq);
    ddc::parallel_fill(kinsp_velocity, kinsp_velocity_eq);
    ddc::parallel_fill(kinsp_temperature, kinsp_temperature_eq);

    neutralsolver.get_derivative(
            derivative,
            neutrals,
            get_const_field(kinsp_density),
            get_const_field(kinsp_velocity),
            get_const_field(kinsp_temperature));

    auto derivative_host = ddc::create_mirror_view_and_copy(derivative);

    double error_l1(0);
    double max_derivative(0);
    ddc::for_each(get_idx_range(neutrals), [&](IdxSpMomX const ispmx) {
        double const neutral_val(neutrals_init_host(ispmx));

        double const lx_2((x_max + x_min) / 2.);
        CoordX coordx(ddc::coordinate(ddc::select<GridX>(ispmx)));
        double const neutral_val_deriv((lx_2 - coordx) * neutral_val);
        double const neutral_val_deriv2(neutral_val * ((lx_2 - coordx) * (lx_2 - coordx) - 1.));

        double const mass_ratio(mass(my_ielec) / mass(my_iion));

        double const advec_term(
                std::sqrt(mass_ratio) * kinsp_velocity_eq * charge_exchange_val
                / (charge_exchange_val + ionization_val));
        double const diffusive_term(
                normalization_coeff * kinsp_temperature_eq
                / (neutral_mass * kinsp_density_eq * (charge_exchange_val + ionization_val)));
        double const derivative_pred(
                -advec_term * neutral_val_deriv + diffusive_term * neutral_val_deriv2);

        error_l1 += abs(derivative_pred - derivative_host(ispmx));
        max_derivative = abs(derivative_pred) > max_derivative ? derivative_pred : max_derivative;
    });

    EXPECT_LE(error_l1 / max_derivative, 0.01);
}
