// SPDX-License-Identifier: MIT
#include <cmath>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <collisions_inter.hpp>
#include <collisions_utils.hpp>
#include <fluid_moments.hpp>
#include <geometry.hpp>
#include <irighthandside.hpp>
#include <maxwellianequilibrium.hpp>
#include <pdi.h>
#include <quadrature.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

TEST(CollisionsInter, CollisionsInter)
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
    double const mass_ion(400.), mass_elec(1.);
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

    std::vector<double> deltat_list = {0.1, 0.01};
    std::vector<double> error_deltat;
    for (double deltat : deltat_list) {
        // Initialization of the distribution function as a maxwellian with a
        // different electron and ion temperatures
        double const density_init(1.);
        DFieldSp temperature_init(dom_sp);
        temperature_init(my_iion) = 1.;
        temperature_init(my_ielec) = 1.2;
        double const fluid_velocity_init(0.);
        ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
            device_t<DFieldVx> finit_device(gridvx);
            MaxwellianEquilibrium::compute_maxwellian(
                    finit_device.span_view(),
                    density_init,
                    temperature_init(ddc::select<IDimSp>(ispx)),
                    fluid_velocity_init);
            auto finit = ddc::create_mirror_view_and_copy(finit_device.span_view());
            ddc::deepcopy(allfdistribu_device[ispx], finit);
        });


        double const nustar0(0.1);
        CollisionsInter collisions(mesh, nustar0);

        DFieldVx const quadrature_coeffs = trapezoid_quadrature_coefficients(gridvx);
        Quadrature<IDimVx> integrate(quadrature_coeffs);
        FluidMoments moments(integrate);

        auto allfdistribu = ddc::create_mirror_view_and_copy(allfdistribu_device.span_view());
        DFieldSpX nustar_profile(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
        compute_nustar_profile(nustar_profile.span_view(), nustar0);

        int const nbiter(100);
        for (int iter(0); iter < nbiter; iter++) {
            collisions(allfdistribu_device, deltat);
        }
        ddc::deepcopy(allfdistribu, allfdistribu_device);

        double error_L1(0);
        ddc::for_each(gridx, [&](IndexX const ix) {
            DFieldSp density(dom_sp);
            DFieldSp fluid_velocity(dom_sp);
            DFieldSp temperature(dom_sp);
            ddc::for_each(dom_sp, [&](IndexSp const isp) {
                IndexSpX const ispx(isp, ix);
                moments(density(isp), allfdistribu[ispx], FluidMoments::s_density);
                moments(fluid_velocity(isp),
                        allfdistribu[ispx],
                        density(isp),
                        FluidMoments::s_velocity);
                moments(temperature(isp),
                        allfdistribu[ispx],
                        density(isp),
                        fluid_velocity(isp),
                        FluidMoments::s_temperature);
            });

            //Collision frequencies, momentum and energy exchange terms
            DFieldSp nustar_profile_copy(nustar_profile[ix]);
            DFieldSp collfreq_ab(dom_sp);
            compute_collfreq_ab(collfreq_ab.span_view(), nustar_profile_copy, density, temperature);

            double const me_on_memi(mass(my_ielec) / (mass(my_ielec) + mass(my_iion)));
            // test : dlog(T_e - T_i)/dt = -12nu_ei*m_e/(m_e+m_b)
            // should be verified
            double const error = std::fabs(
                    std::log(std::fabs(temperature(my_ielec) - temperature(my_iion)))
                    - std::log(std::fabs(temperature_init(my_ielec) - temperature_init(my_iion)))
                    + 12 * collfreq_ab(my_ielec) * me_on_memi * nbiter * deltat);
            error_L1 += error;
        });
        error_L1 = error_L1 / x_size;
        error_deltat.emplace_back(error_L1);
    }

    int const order = std::log(deltat_list.front() / deltat_list.back());
    double relative_error = std::fabs(
            (error_deltat.back() * pow(10., order) - error_deltat.front()) / error_deltat.front());
    // relative error should be less than 5%
    EXPECT_LE(relative_error, 0.05);

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
