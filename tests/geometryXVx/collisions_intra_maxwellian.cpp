// SPDX-License-Identifier: MIT
#define _USE_MATH_DEFINES

#include <cmath>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <collisions_intra.hpp>
#include <collisions_utils.hpp>
#include <fluid_moments.hpp>
#include <geometry.hpp>
#include <irighthandside.hpp>
#include <maxwellianequilibrium.hpp>
#include <pdi.h>
#include <quadrature.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

/**
 * Intra species collisions applied on a maxwellian should not change the distribution function
 */
TEST(CollisionsIntraMaxwellian, CollisionsIntraMaxwellian)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_size(5);

    CoordVx const vx_min(-10);
    CoordVx const vx_max(10);
    IVectVx const vx_size(600);

    IVectSp const nb_kinspecies(2);

    IDomainSp const dom_sp(IndexSp(0), nb_kinspecies);
    IndexSp const my_iion = dom_sp.front();
    IndexSp const my_ielec = dom_sp.back();

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());

    IDomainX interpolation_domain_x(SplineInterpPointsX::get_domain());
    IDomainVx interpolation_domain_vx(SplineInterpPointsVx::get_domain());

    SplineXBuilder const builder_x(interpolation_domain_x);

    SplineVxBuilder const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSpXVx const mesh(dom_sp, gridx, gridvx);

    FieldSp<int> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    DFieldSp masses(dom_sp);
    double const mass_ion(400), mass_elec(1);
    masses(my_ielec) = mass_elec;
    masses(my_iion) = mass_ion;
    FieldSp<int> init_perturb_mode(dom_sp);
    ddc::fill(init_perturb_mode, 0);
    DFieldSp init_perturb_amplitude(dom_sp);
    ddc::fill(init_perturb_amplitude, 0);

    // Initialization of the distribution function as a maxwellian
    ddc::init_discrete_space<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));
    device_t<DFieldSpXVx> allfdistribu_device(mesh);

    // Initialization of the distribution function as a maxwellian with
    // moments depending on space
    DFieldSpX density_init(ddc::select<IDimSp, IDimX>(mesh));
    DFieldSpX mean_velocity_init(ddc::select<IDimSp, IDimX>(mesh));
    DFieldSpX temperature_init(ddc::select<IDimSp, IDimX>(mesh));
    ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
        double const density = 1.;
        double const density_ampl = 0.1;
        double const mean_velocity = 0.;
        double const mean_velocity_ampl = 0.2;
        double const temperature = 1;
        double const temperature_ampl = 0.3;

        double const coordx = ddc::coordinate(ddc::select<IDimX>(ispx));
        density_init(ispx)
                = density
                  + density_ampl * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        mean_velocity_init(ispx)
                = mean_velocity
                  + mean_velocity_ampl
                            * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        temperature_init(ispx)
                = temperature
                  + temperature_ampl * std::sin(2 * M_PI * coordx / ddc::coordinate(gridx.back()));
        device_t<DFieldVx> finit_device(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                finit_device.span_view(),
                density_init(ispx),
                temperature_init(ispx),
                mean_velocity_init(ispx));
        auto finit = ddc::create_mirror_view_and_copy(finit_device.span_view());
        ddc::deepcopy(allfdistribu_device[ispx], finit);
    });
    auto allfdistribu = ddc::create_mirror_view_and_copy(allfdistribu_device.span_view());

    double const nustar0(0.1);
    double const deltat(0.1);
    CollisionsIntra collisions(mesh, nustar0);

    // test of the get_elec_index
    EXPECT_EQ(charge(ielec()), -1);

    // nustar profile
    DFieldSpX nustar_profile(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    compute_nustar_profile(nustar_profile.span_view(), nustar0);
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        if (charge(ddc::select<IDimSp>(ispx)) < 0.) {
            double const pred(1 / x_max * nustar0);
            EXPECT_LE(std::fabs(nustar_profile(ispx) - pred), 1e-12);
        } else {
            double const pred(1 / (std::sqrt(mass_ion) * x_max) * nustar0);
            EXPECT_LE(std::fabs(nustar_profile(ispx) - pred), 1e-12);
        }
    });

    //collfreq
    DFieldSpX collfreq(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    compute_collfreq(
            collfreq.span_view(),
            nustar_profile.span_cview(),
            density_init.span_cview(),
            temperature_init.span_cview());

    ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
        if (charge(ddc::select<IDimSp>(ispx)) < 0.) {
            double const pred(
                    1 / x_max * nustar0 * density_init(ispx)
                    / std::pow(temperature_init(ispx), 1.5));
            EXPECT_LE(std::fabs(collfreq(ispx) - pred), 1e-12);
        } else {
            double const pred(
                    1 / (std::sqrt(mass_ion) * x_max) * nustar0 * density_init(ispx)
                    / std::pow(temperature_init(ispx), 1.5));
            EXPECT_LE(std::fabs(collfreq(ispx) - pred), 1e-12);
        }
    });

    // diffusion coefficient
    ddc::Chunk<
            double,
            ddc::DiscreteDomain<IDimSp, IDimX, CollisionsIntra::ghosted_vx_point_sampling>>
            Dcoll(collisions.get_mesh_ghosted());
    compute_Dcoll<CollisionsIntra::ghosted_vx_point_sampling>(
            Dcoll.span_view(),
            collfreq.span_cview(),
            density_init.span_cview(),
            temperature_init.span_cview());

    ddc::Chunk<
            double,
            ddc::DiscreteDomain<IDimSp, IDimX, CollisionsIntra::ghosted_vx_point_sampling>>
            dvDcoll(collisions.get_mesh_ghosted());
    compute_dvDcoll<CollisionsIntra::ghosted_vx_point_sampling>(
            dvDcoll.span_view(),
            collfreq.span_cview(),
            density_init.span_cview(),
            temperature_init.span_cview());

    // kernel maxwellian fluid moments
    DFieldSpX Vcoll(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX Tcoll(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    compute_Vcoll_Tcoll(
            Vcoll.span_view(),
            Tcoll.span_view(),
            allfdistribu.span_cview(),
            Dcoll.span_cview(),
            dvDcoll.span_cview());

    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        EXPECT_LE(std::fabs(Vcoll(ispx) - mean_velocity_init(ispx)), 1e-12);
        EXPECT_LE(std::fabs(Tcoll(ispx) - temperature_init(ispx)), 1e-12);
    });

    collisions(allfdistribu_device, deltat);
    ddc::deepcopy(allfdistribu, allfdistribu_device);

    // collision operator should not change densiy, mean_velocity and temperature
    DFieldSpX density_res(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX mean_velocity_res(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX temperature_res(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldVx const quadrature_coeffs
            = trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu));
    Quadrature<IDimVx> integrate(quadrature_coeffs);
    FluidMoments moments(integrate);

    moments(density_res.span_view(), allfdistribu.span_cview(), FluidMoments::s_density);
    moments(mean_velocity_res.span_view(),
            allfdistribu.span_cview(),
            density_res.span_cview(),
            FluidMoments::s_velocity);
    moments(temperature_res.span_view(),
            allfdistribu.span_cview(),
            density_res.span_cview(),
            mean_velocity_res.span_cview(),
            FluidMoments::s_temperature);

    double const tol = 1.e-6;
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        EXPECT_LE(std::fabs(density_res(ispx) - density_init(ispx)), tol);
        EXPECT_LE(std::fabs(mean_velocity_res(ispx) - mean_velocity_init(ispx)), tol);
        EXPECT_LE(std::fabs(temperature_res(ispx) - temperature_init(ispx)), tol);
    });

    // * Intra species collisions applied on a perturbed distribution function
    // * should make it relax towards a maxwellian, i.e. Vcoll = mean_velocity, Tcoll = temperature at equilibrium
    ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
        double const density = 1.;
        double const density_ampl = 0.1;
        double const mean_velocity = 0.;
        double const mean_velocity_ampl = 0.2;
        double const temperature = 1;
        double const temperature_ampl = 0.3;
        double const coordx = ddc::coordinate(ddc::select<IDimX>(ispx));
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
        device_t<DFieldVx> finit_device(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                finit_device.span_view(),
                density_loc,
                temperature_loc,
                mean_velocity_loc);
        auto finit = ddc::create_mirror_view_and_copy(finit_device.span_view());
        ddc::deepcopy(allfdistribu[ispx], finit);
    });

    // initial perturbation
    double const perturb_velocity = 0.5;
    double const perturb_temperature = 0.3;
    double const perturb_density = 0.5;
    double const coeff_mxw(perturb_density / std::sqrt(2 * M_PI * perturb_temperature));

    ddc::for_each(mesh, [&](IndexSpXVx const ispxvx) {
        CoordVx const coordvx = ddc::coordinate(ddc::select<IDimVx>(ispxvx));
        double const coordvx_sq = (coordvx - perturb_velocity) * (coordvx - perturb_velocity);
        double const perturb
                = coeff_mxw * std::exp(-coordvx_sq * coordvx_sq / (2 * perturb_temperature));

        allfdistribu(ispxvx) = allfdistribu(ispxvx) + perturb;
    });
    ddc::deepcopy(allfdistribu_device, allfdistribu);

    // before the collisions the perturbed distribution should have T != Tcoll and V != Vcoll
    // this test is performed on electrons since the largest error is made for them
    // (lightest species constraining the timestep)
    moments(density_res.span_view(), allfdistribu.span_cview(), FluidMoments::s_density);
    moments(mean_velocity_res.span_view(),
            allfdistribu.span_cview(),
            density_res.span_cview(),
            FluidMoments::s_velocity);
    moments(temperature_res.span_view(),
            allfdistribu.span_cview(),
            density_res.span_cview(),
            mean_velocity_res.span_cview(),
            FluidMoments::s_temperature);


    ddc::for_each(ddc::get_domain<IDimX>(allfdistribu), [&](IndexX const ix) {
        EXPECT_GE(std::fabs(mean_velocity_res(ielec(), ix) - Vcoll(ielec(), ix)), 1.e-4);
        EXPECT_GE(std::fabs(temperature_res(ielec(), ix) - Tcoll(ielec(), ix)), 1.e-4);
    });


    int const nbsteps = 300;
    for (int iter = 0; iter < nbsteps; ++iter) {
        collisions(allfdistribu_device, deltat);
    };
    ddc::deepcopy(allfdistribu, allfdistribu_device);

    // Vcoll and Tcoll calculation
    compute_collfreq(
            collfreq.span_view(),
            nustar_profile.span_cview(),
            density_init.span_cview(),
            temperature_init.span_cview());

    compute_Dcoll<CollisionsIntra::ghosted_vx_point_sampling>(
            Dcoll.span_view(),
            collfreq.span_cview(),
            density_init.span_cview(),
            temperature_init.span_cview());

    compute_dvDcoll<CollisionsIntra::ghosted_vx_point_sampling>(
            dvDcoll.span_view(),
            collfreq.span_cview(),
            density_init.span_cview(),
            temperature_init.span_cview());

    compute_Vcoll_Tcoll(
            Vcoll.span_view(),
            Tcoll.span_view(),
            allfdistribu.span_cview(),
            Dcoll.span_cview(),
            dvDcoll.span_cview());

    moments(density_res.span_view(), allfdistribu.span_cview(), FluidMoments::s_density);
    moments(mean_velocity_res.span_view(),
            allfdistribu.span_cview(),
            density_res.span_cview(),
            FluidMoments::s_velocity);
    moments(temperature_res.span_view(),
            allfdistribu.span_cview(),
            density_res.span_cview(),
            mean_velocity_res.span_cview(),
            FluidMoments::s_temperature);

    ddc::for_each(ddc::get_domain<IDimX>(allfdistribu), [&](IndexX const ix) {
        EXPECT_LE(std::fabs(mean_velocity_res(ielec(), ix) - Vcoll(ielec(), ix)), 1.e-4);
        EXPECT_LE(std::fabs(temperature_res(ielec(), ix) - Tcoll(ielec(), ix)), 1.e-4);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
