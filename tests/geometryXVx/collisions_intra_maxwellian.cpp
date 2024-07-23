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

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());

    IDomainX interpolation_domain_x(SplineInterpPointsX::get_domain<IDimX>());
    IDomainVx interpolation_domain_vx(SplineInterpPointsVx::get_domain<IDimVx>());

    SplineXBuilder_1d const builder_x(interpolation_domain_x);

    SplineVxBuilder_1d const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSpXVx const mesh(dom_sp, gridx, gridvx);

    host_t<DFieldSp> charges(dom_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldSp> masses(dom_sp);
    double const mass_ion(400.), mass_elec(1.);
    masses(my_ielec) = mass_elec;
    masses(my_iion) = mass_ion;

    // Initialization of the distribution function as a maxwellian
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));
    DFieldSpXVx allfdistribu(mesh);

    // Initialization of the distribution function as a maxwellian with
    // moments depending on space
    host_t<DFieldSpX> density_init_host(ddc::select<IDimSp, IDimX>(mesh));
    host_t<DFieldSpX> mean_velocity_init_host(ddc::select<IDimSp, IDimX>(mesh));
    host_t<DFieldSpX> temperature_init_host(ddc::select<IDimSp, IDimX>(mesh));
    ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
        double const density = 1.;
        double const density_ampl = 0.1;
        double const mean_velocity = 0.;
        double const mean_velocity_ampl = 0.2;
        double const temperature = 1;
        double const temperature_ampl = 0.3;

        double const coordx = ddc::coordinate(ddc::select<IDimX>(ispx));
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
        DFieldVx finit(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                finit.span_view(),
                density_init_host(ispx),
                temperature_init_host(ispx),
                mean_velocity_init_host(ispx));
        auto finit_host = ddc::create_mirror_view_and_copy(finit.span_view());
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu.span_view());

    double const nustar0(0.1);
    double const deltat(0.1);
    CollisionsIntra collisions(mesh, nustar0);

    // test of the get_elec_index
    EXPECT_EQ(charge(ielec()), -1.);

    // nustar profile
    DFieldSpX nustar_profile_alloc(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    DSpanSpX nustar_profile(nustar_profile_alloc.span_view());
    compute_nustar_profile(nustar_profile, nustar0);

    host_t<DFieldSpX> nustar_profile_host(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    ddc::parallel_deepcopy(nustar_profile_host, nustar_profile);
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host), [&](IndexSpX const ispx) {
        if (charge(ddc::select<IDimSp>(ispx)) < 0.) {
            double const pred(1 / x_max * nustar0);
            EXPECT_LE(std::fabs(nustar_profile_host(ispx) - pred), 1e-12);
        } else {
            double const pred(1 / (std::sqrt(mass_ion) * x_max) * nustar0);
            EXPECT_LE(std::fabs(nustar_profile_host(ispx) - pred), 1e-12);
        }
    });

    //collfreq
    DFieldSpX collfreq_f(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    auto collfreq = collfreq_f.span_view();

    DFieldSpX density_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    ddc::parallel_deepcopy(density_init, density_init_host);

    DFieldSpX temperature_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    ddc::parallel_deepcopy(temperature_init, temperature_init_host);

    compute_collfreq(collfreq, nustar_profile, density_init, temperature_init);

    host_t<DFieldSpX> collfreq_host(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    ddc::parallel_deepcopy(collfreq_host, collfreq);

    ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
        if (charge(ddc::select<IDimSp>(ispx)) < 0.) {
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
    device_t<ddc::Chunk<double, ddc::DiscreteDomain<IDimSp, IDimX, CollisionsIntra::GhostedVx>>>
            Dcoll_f(collisions.get_mesh_ghosted());
    auto Dcoll = Dcoll_f.span_view();
    compute_Dcoll<CollisionsIntra::GhostedVx>(Dcoll, collfreq, density_init, temperature_init);

    device_t<ddc::Chunk<double, ddc::DiscreteDomain<IDimSp, IDimX, CollisionsIntra::GhostedVx>>>
            dvDcoll_f(collisions.get_mesh_ghosted());
    auto dvDcoll = dvDcoll_f.span_view();
    compute_dvDcoll<CollisionsIntra::GhostedVx>(dvDcoll, collfreq, density_init, temperature_init);

    // kernel maxwellian fluid moments
    DFieldSpX Vcoll_f(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    DFieldSpX Tcoll_f(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    auto Vcoll = Vcoll_f.span_view();
    auto Tcoll = Tcoll_f.span_view();
    compute_Vcoll_Tcoll<CollisionsIntra::GhostedVx>(Vcoll, Tcoll, allfdistribu, Dcoll, dvDcoll);

    host_t<DFieldSpX> Vcoll_host(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    host_t<DFieldSpX> Tcoll_host(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    ddc::parallel_deepcopy(Vcoll_host, Vcoll);
    ddc::parallel_deepcopy(Tcoll_host, Tcoll);

    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host), [&](IndexSpX const ispx) {
        EXPECT_LE(std::fabs(Vcoll_host(ispx) - mean_velocity_init_host(ispx)), 1e-12);
        EXPECT_LE(std::fabs(Tcoll_host(ispx) - temperature_init_host(ispx)), 1e-12);
    });

    collisions(allfdistribu, deltat);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    // collision operator should not change densiy, mean_velocity and temperature
    DFieldSpX density_res(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    DFieldSpX mean_velocity_res(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    DFieldSpX temperature_res(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));

    DFieldVx const quadrature_coeffs
            = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridvx);
    Quadrature<IDomainVx, IDomainSpXVx> integrate(quadrature_coeffs.span_cview());

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
    auto mean_velocity_res_host = ddc::create_mirror_view_and_copy(mean_velocity_res.span_view());
    auto temperature_res_host = ddc::create_mirror_view_and_copy(temperature_res.span_view());
    auto density_res_host = ddc::create_mirror_view_and_copy(density_res.span_view());

    double const tol = 1.e-6;
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host), [&](IndexSpX const ispx) {
        EXPECT_LE(std::fabs(density_res_host(ispx) - density_init_host(ispx)), tol);
        EXPECT_LE(std::fabs(mean_velocity_res_host(ispx) - mean_velocity_init_host(ispx)), tol);
        EXPECT_LE(std::fabs(temperature_res_host(ispx) - temperature_init_host(ispx)), tol);
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
        DFieldVx finit(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                finit.span_view(),
                density_loc,
                temperature_loc,
                mean_velocity_loc);
        auto finit_host = ddc::create_mirror_view_and_copy(finit.span_view());
        ddc::parallel_deepcopy(allfdistribu_host[ispx], finit_host);
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

        allfdistribu_host(ispxvx) = allfdistribu_host(ispxvx) + perturb;
    });
    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);

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
    ddc::parallel_deepcopy(mean_velocity_res_host, mean_velocity_res);
    ddc::parallel_deepcopy(temperature_res_host, temperature_res);

    ddc::for_each(ddc::get_domain<IDimX>(allfdistribu_host), [&](IndexX const ix) {
        EXPECT_GE(std::fabs(mean_velocity_res_host(ielec(), ix) - Vcoll_host(ielec(), ix)), 1.e-4);
        EXPECT_GE(std::fabs(temperature_res_host(ielec(), ix) - Tcoll_host(ielec(), ix)), 1.e-4);
    });

    int const nbsteps = 300;
    for (int iter = 0; iter < nbsteps; ++iter) {
        collisions(allfdistribu, deltat);
    };
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    // Vcoll and Tcoll calculation
    compute_collfreq(collfreq, nustar_profile, density_init, temperature_init);

    compute_Dcoll<CollisionsIntra::GhostedVx>(Dcoll, collfreq, density_init, temperature_init);

    compute_dvDcoll<CollisionsIntra::GhostedVx>(dvDcoll, collfreq, density_init, temperature_init);

    compute_Vcoll_Tcoll<CollisionsIntra::GhostedVx>(Vcoll, Tcoll, allfdistribu, Dcoll, dvDcoll);

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

    ddc::parallel_deepcopy(Vcoll_host, Vcoll);
    ddc::parallel_deepcopy(Tcoll_host, Tcoll);
    ddc::parallel_deepcopy(mean_velocity_res_host, mean_velocity_res);
    ddc::parallel_deepcopy(temperature_res_host, temperature_res);
    ddc::for_each(ddc::get_domain<IDimX>(allfdistribu_host), [&](IndexX const ix) {
        EXPECT_LE(std::fabs(mean_velocity_res_host(ielec(), ix) - Vcoll_host(ielec(), ix)), 1.e-4);
        EXPECT_LE(std::fabs(temperature_res_host(ielec(), ix) - Tcoll_host(ielec(), ix)), 1.e-4);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
