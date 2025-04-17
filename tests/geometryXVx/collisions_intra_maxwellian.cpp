// SPDX-License-Identifier: MIT
#include "ddc_alias_inline_functions.hpp"
#define _USE_MATH_DEFINES

#include <cmath>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>

#include "collisions_intra.hpp"
#include "collisions_utils.hpp"
#include "fluid_moments.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "maxwellianequilibrium.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"
#include "trapezoid_quadrature.hpp"

/**
 * Intra species collisions applied on a maxwellian should not change the distribution function
 */
TEST(CollisionsIntraMaxwellian, CollisionsIntraMaxwellian)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(5);

    CoordVx const vx_min(-10);
    CoordVx const vx_max(10);
    IdxStepVx const vx_size(600);

    IdxStepSp const nb_kinspecies(2);

    IdxRangeSp const idx_range_sp(IdxSp(0), nb_kinspecies);
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

    IdxRangeSpXVx const mesh(idx_range_sp, gridx, gridvx);

    host_t<DFieldMemSp> charges(idx_range_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldMemSp> masses(idx_range_sp);
    double const mass_ion(400.), mass_elec(1.);
    masses(my_ielec) = mass_elec;
    masses(my_iion) = mass_ion;

    // Initialisation of the distribution function as a maxwellian
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
    DFieldMemSpXVx allfdistribu(mesh);

    // Initialisation of the distribution function as a maxwellian with
    // moments depending on space
    host_t<DFieldMemSpX> density_init_host(ddc::select<Species, GridX>(mesh));
    host_t<DFieldMemSpX> mean_velocity_init_host(ddc::select<Species, GridX>(mesh));
    host_t<DFieldMemSpX> temperature_init_host(ddc::select<Species, GridX>(mesh));
    ddc::for_each(ddc::select<Species, GridX>(mesh), [&](IdxSpX const ispx) {
        double const density = 1.;
        double const density_ampl = 0.1;
        double const mean_velocity = 0.;
        double const mean_velocity_ampl = 0.2;
        double const temperature = 1;
        double const temperature_ampl = 0.3;

        double const coordx = ddc::coordinate(ddc::select<GridX>(ispx));
        density_init_host(ispx)
                = density
                  + density_ampl * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        mean_velocity_init_host(ispx)
                = mean_velocity
                  + mean_velocity_ampl
                            * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        temperature_init_host(ispx)
                = temperature
                  + temperature_ampl * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        DFieldMemVx finit(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                get_field(finit),
                density_init_host(ispx),
                temperature_init_host(ispx),
                mean_velocity_init_host(ispx));
        auto finit_host = ddc::create_mirror_view_and_copy(get_field(finit));
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));

    double const nustar0(0.1);
    double const deltat(0.1);
    CollisionsIntra collisions(mesh, nustar0);

    // test of the get_elec_index
    EXPECT_EQ(charge(ielec()), -1.);

    // nustar profile
    DFieldMemSpX nustar_profile_alloc(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldSpX nustar_profile(get_field(nustar_profile_alloc));
    compute_nustar_profile(nustar_profile, nustar0);

    host_t<DFieldMemSpX> nustar_profile_host(get_idx_range<Species, GridX>(allfdistribu_host));
    ddc::parallel_deepcopy(nustar_profile_host, nustar_profile);
    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu_host), [&](IdxSpX const ispx) {
        if (charge(ddc::select<Species>(ispx)) < 0.) {
            double const pred(1 / x_max * nustar0);
            EXPECT_LE(std::fabs(nustar_profile_host(ispx) - pred), 1e-12);
        } else {
            double const pred(1 / (std::sqrt(mass_ion) * x_max) * nustar0);
            EXPECT_LE(std::fabs(nustar_profile_host(ispx) - pred), 1e-12);
        }
    });

    //collfreq
    DFieldMemSpX collfreq_alloc(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldSpX collfreq = get_field(collfreq_alloc);

    DFieldMemSpX density_init(get_idx_range<Species, GridX>(allfdistribu_host));
    ddc::parallel_deepcopy(density_init, density_init_host);

    DFieldMemSpX temperature_init(get_idx_range<Species, GridX>(allfdistribu_host));
    ddc::parallel_deepcopy(temperature_init, temperature_init_host);

    compute_collfreq(
            collfreq,
            get_const_field(nustar_profile),
            get_const_field(density_init),
            get_const_field(temperature_init));

    host_t<DFieldMemSpX> collfreq_host(get_idx_range<Species, GridX>(allfdistribu_host));
    ddc::parallel_deepcopy(collfreq_host, collfreq);

    ddc::for_each(ddc::select<Species, GridX>(mesh), [&](IdxSpX const ispx) {
        if (charge(ddc::select<Species>(ispx)) < 0.) {
            double const pred(
                    1 / x_max * nustar0 * density_init_host(ispx)
                    / std::pow(temperature_init_host(ispx), 1.5));
            EXPECT_LE(std::fabs(collfreq_host(ispx) - pred), 1e-12);
        } else {
            double const pred(
                    1 / (std::sqrt(mass_ion) * x_max) * nustar0 * density_init_host(ispx)
                    / std::pow(temperature_init_host(ispx), 1.5));
            EXPECT_LE(std::fabs(collfreq_host(ispx) - pred), 1e-12);
        }
    });

    // diffusion coefficient
    DFieldMem<IdxRange<Species, GridX, CollisionsIntra::GhostedVx>> Dcoll_alloc(
            collisions.get_mesh_ghosted());
    DField<IdxRange<Species, GridX, CollisionsIntra::GhostedVx>> Dcoll = get_field(Dcoll_alloc);
    compute_Dcoll<CollisionsIntra::GhostedVx>(
            Dcoll,
            get_const_field(collfreq),
            get_const_field(density_init),
            get_const_field(temperature_init));

    DFieldMem<IdxRange<Species, GridX, CollisionsIntra::GhostedVx>> dvDcoll_alloc(
            collisions.get_mesh_ghosted());
    DField<IdxRange<Species, GridX, CollisionsIntra::GhostedVx>> dvDcoll = get_field(dvDcoll_alloc);
    compute_dvDcoll<CollisionsIntra::GhostedVx>(
            dvDcoll,
            get_const_field(collfreq),
            get_const_field(density_init),
            get_const_field(temperature_init));

    // kernel maxwellian fluid moments
    DFieldMemSpX Vcoll_alloc(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldMemSpX Tcoll_alloc(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldSpX Vcoll = get_field(Vcoll_alloc);
    DFieldSpX Tcoll = get_field(Tcoll_alloc);
    compute_Vcoll_Tcoll<CollisionsIntra::GhostedVx>(
            Vcoll,
            get_field(Tcoll),
            get_const_field(allfdistribu),
            get_field(Dcoll),
            get_field(dvDcoll));

    host_t<DFieldMemSpX> Vcoll_host(get_idx_range<Species, GridX>(allfdistribu_host));
    host_t<DFieldMemSpX> Tcoll_host(get_idx_range<Species, GridX>(allfdistribu_host));
    ddc::parallel_deepcopy(Vcoll_host, Vcoll);
    ddc::parallel_deepcopy(Tcoll_host, Tcoll);

    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu_host), [&](IdxSpX const ispx) {
        EXPECT_LE(std::fabs(Vcoll_host(ispx) - mean_velocity_init_host(ispx)), 1e-12);
        EXPECT_LE(std::fabs(Tcoll_host(ispx) - temperature_init_host(ispx)), 1e-12);
    });

    collisions(get_field(allfdistribu), deltat);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    // collision operator should not change densiy, mean_velocity and temperature
    DFieldMemSpX density_res(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldMemSpX mean_velocity_res(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldMemSpX temperature_res(get_idx_range<Species, GridX>(allfdistribu_host));

    DFieldMemVx const quadrature_coeffs
            = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridvx);
    Quadrature<IdxRangeVx, IdxRangeSpXVx> integrate(get_const_field(quadrature_coeffs));

    FluidMoments moments(integrate);

    moments(get_field(density_res), get_const_field(allfdistribu), FluidMoments::s_density);
    moments(get_field(mean_velocity_res),
            get_const_field(allfdistribu),
            get_const_field(density_res),
            FluidMoments::s_velocity);
    moments(get_field(temperature_res),
            get_const_field(allfdistribu),
            get_const_field(density_res),
            get_const_field(mean_velocity_res),
            FluidMoments::s_temperature);
    auto mean_velocity_res_host = ddc::create_mirror_view_and_copy(get_field(mean_velocity_res));
    auto temperature_res_host = ddc::create_mirror_view_and_copy(get_field(temperature_res));
    auto density_res_host = ddc::create_mirror_view_and_copy(get_field(density_res));

    double const tol = 1.e-6;
    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu_host), [&](IdxSpX const ispx) {
        EXPECT_LE(std::fabs(density_res_host(ispx) - density_init_host(ispx)), tol);
        EXPECT_LE(std::fabs(mean_velocity_res_host(ispx) - mean_velocity_init_host(ispx)), tol);
        EXPECT_LE(std::fabs(temperature_res_host(ispx) - temperature_init_host(ispx)), tol);
    });

    // * Intra species collisions applied on a perturbed distribution function
    // * should make it relax towards a maxwellian, i.e. Vcoll = mean_velocity, Tcoll = temperature at equilibrium
    ddc::for_each(ddc::select<Species, GridX>(mesh), [&](IdxSpX const ispx) {
        double const density = 1.;
        double const density_ampl = 0.1;
        double const mean_velocity = 0.;
        double const mean_velocity_ampl = 0.2;
        double const temperature = 1;
        double const temperature_ampl = 0.3;
        double const coordx = ddc::coordinate(ddc::select<GridX>(ispx));
        double const density_loc
                = density
                  + density_ampl * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        double const mean_velocity_loc
                = mean_velocity
                  + mean_velocity_ampl
                            * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        double const temperature_loc
                = temperature
                  + temperature_ampl * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        DFieldMemVx finit(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                get_field(finit),
                density_loc,
                temperature_loc,
                mean_velocity_loc);
        auto finit_host = ddc::create_mirror_view_and_copy(get_field(finit));
        ddc::parallel_deepcopy(allfdistribu_host[ispx], finit_host);
    });

    // initial perturbation
    double const perturb_velocity = 0.5;
    double const perturb_temperature = 0.3;
    double const perturb_density = 0.5;
    double const coeff_mxw(perturb_density / std::sqrt(2 * M_PI * perturb_temperature));

    ddc::for_each(mesh, [&](IdxSpXVx const ispxvx) {
        CoordVx const coordvx = ddc::coordinate(ddc::select<GridVx>(ispxvx));
        double const coordvx_sq = (coordvx - perturb_velocity) * (coordvx - perturb_velocity);
        double const perturb
                = coeff_mxw * std::exp(-coordvx_sq * coordvx_sq / (2 * perturb_temperature));

        allfdistribu_host(ispxvx) = allfdistribu_host(ispxvx) + perturb;
    });
    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);

    // before the collisions the perturbed distribution should have T != Tcoll and V != Vcoll
    // this test is performed on electrons since the largest error is made for them
    // (lightest species constraining the timestep)
    moments(get_field(density_res), get_const_field(allfdistribu), FluidMoments::s_density);
    moments(get_field(mean_velocity_res),
            get_const_field(allfdistribu),
            get_const_field(density_res),
            FluidMoments::s_velocity);
    moments(get_field(temperature_res),
            get_const_field(allfdistribu),
            get_const_field(density_res),
            get_const_field(mean_velocity_res),
            FluidMoments::s_temperature);
    ddc::parallel_deepcopy(mean_velocity_res_host, mean_velocity_res);
    ddc::parallel_deepcopy(temperature_res_host, temperature_res);

    ddc::for_each(get_idx_range<GridX>(allfdistribu_host), [&](IdxX const ix) {
        EXPECT_GE(std::fabs(mean_velocity_res_host(ielec(), ix) - Vcoll_host(ielec(), ix)), 1.e-4);
        EXPECT_GE(std::fabs(temperature_res_host(ielec(), ix) - Tcoll_host(ielec(), ix)), 1.e-4);
    });

    int const nbsteps = 300;
    for (int iter = 0; iter < nbsteps; ++iter) {
        collisions(get_field(allfdistribu), deltat);
    };
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    // Vcoll and Tcoll calculation
    compute_collfreq(
            collfreq,
            get_const_field(nustar_profile),
            get_const_field(density_init),
            get_const_field(temperature_init));

    compute_Dcoll<CollisionsIntra::GhostedVx>(
            Dcoll,
            get_const_field(collfreq),
            get_const_field(density_init),
            get_const_field(temperature_init));

    compute_dvDcoll<CollisionsIntra::GhostedVx>(
            dvDcoll,
            get_const_field(collfreq),
            get_const_field(density_init),
            get_const_field(temperature_init));

    compute_Vcoll_Tcoll<CollisionsIntra::GhostedVx>(
            Vcoll,
            get_field(Tcoll),
            get_const_field(allfdistribu),
            get_field(Dcoll),
            get_field(dvDcoll));

    moments(get_field(density_res), get_const_field(allfdistribu), FluidMoments::s_density);
    moments(get_field(mean_velocity_res),
            get_const_field(allfdistribu),
            get_const_field(density_res),
            FluidMoments::s_velocity);
    moments(get_field(temperature_res),
            get_const_field(allfdistribu),
            get_const_field(density_res),
            get_const_field(mean_velocity_res),
            FluidMoments::s_temperature);

    ddc::parallel_deepcopy(Vcoll_host, Vcoll);
    ddc::parallel_deepcopy(Tcoll_host, Tcoll);
    ddc::parallel_deepcopy(mean_velocity_res_host, mean_velocity_res);
    ddc::parallel_deepcopy(temperature_res_host, temperature_res);
    ddc::for_each(get_idx_range<GridX>(allfdistribu_host), [&](IdxX const ix) {
        EXPECT_LE(std::fabs(mean_velocity_res_host(ielec(), ix) - Vcoll_host(ielec(), ix)), 1.e-4);
        EXPECT_LE(std::fabs(temperature_res_host(ielec(), ix) - Tcoll_host(ielec(), ix)), 1.e-4);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
