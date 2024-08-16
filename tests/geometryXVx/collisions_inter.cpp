// SPDX-License-Identifier: MIT
#include <cmath>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>

#include "collisions_inter.hpp"
#include "collisions_utils.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "fluid_moments.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "maxwellianequilibrium.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"
#include "trapezoid_quadrature.hpp"

TEST(CollisionsInter, CollisionsInter)
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

    SplineXBuilder_1d const builder_x(gridx);
    SplineVxBuilder_1d const builder_vx(gridvx);

    IdxRangeSpXVx const mesh(idx_range_sp, gridx, gridvx);

    host_t<DFieldMemSp> charges(idx_range_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldMemSp> masses(idx_range_sp);
    double const mass_ion(400.), mass_elec(1.);
    masses(my_ielec) = mass_elec;
    masses(my_iion) = mass_ion;

    // Initialization of the distribution function as a maxwellian
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
    DFieldMemSpXVx allfdistribu(mesh);

    std::vector<double> deltat_list = {0.1, 0.01};
    std::vector<double> error_deltat;
    for (double deltat : deltat_list) {
        // Initialization of the distribution function as a maxwellian with a
        // different electron and ion temperatures
        double const density_init(1.);
        host_t<DFieldMemSp> temperature_init(idx_range_sp);
        temperature_init(my_iion) = 1.;
        temperature_init(my_ielec) = 1.2;
        double const fluid_velocity_init(0.);
        ddc::for_each(ddc::select<Species, GridX>(mesh), [&](IdxSpX const ispx) {
            DFieldMemVx finit(gridvx);
            MaxwellianEquilibrium::compute_maxwellian(
                    get_field(finit),
                    density_init,
                    temperature_init(ddc::select<Species>(ispx)),
                    fluid_velocity_init);
            auto finit_host = ddc::create_mirror_view_and_copy(get_field(finit));
            ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
        });


        double const nustar0(0.1);
        CollisionsInter collisions(mesh, nustar0);

        DFieldMemVx quadrature_coeffs
                = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridvx);

        Quadrature<IdxRangeVx, IdxRangeSpXVx> integrate(get_const_field(quadrature_coeffs));
        FluidMoments moments(integrate);

        auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));
        DFieldMemSpX nustar_profile_alloc(get_idx_range<Species, GridX>(allfdistribu_host));
        DFieldSpX nustar_profile = get_field(nustar_profile_alloc);
        compute_nustar_profile(nustar_profile, nustar0);

        int const nbiter(100);
        for (int iter(0); iter < nbiter; iter++) {
            collisions(allfdistribu, deltat);
        }
        ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

        double error_L1(0);
        DFieldMemSpX density(get_idx_range<Species, GridX>(allfdistribu_host));
        DFieldMemSpX temperature(get_idx_range<Species, GridX>(allfdistribu_host));
        host_t<DFieldMemSpX> density_host(get_idx_range<Species, GridX>(allfdistribu_host));
        host_t<DFieldMemSpX> fluid_velocity_host(get_idx_range<Species, GridX>(allfdistribu_host));
        host_t<DFieldMemSpX> temperature_host(get_idx_range<Species, GridX>(allfdistribu_host));

        ddc::for_each(get_idx_range<Species, GridX>(allfdistribu_host), [&](IdxSpX const ispx) {
            moments(density_host(ispx), allfdistribu[ispx], FluidMoments::s_density);
            moments(fluid_velocity_host(ispx),
                    allfdistribu[ispx],
                    density_host(ispx),
                    FluidMoments::s_velocity);
            moments(temperature_host(ispx),
                    allfdistribu[ispx],
                    density_host(ispx),
                    fluid_velocity_host(ispx),
                    FluidMoments::s_temperature);
        });
        ddc::parallel_deepcopy(temperature, temperature_host);
        ddc::parallel_deepcopy(density, density_host);

        //Collision frequencies, momentum and energy exchange terms
        DFieldMemSpX collfreq_ab(get_idx_range<Species, GridX>(allfdistribu_host));
        compute_collfreq_ab(get_field(collfreq_ab), nustar_profile, density, temperature);
        auto collfreq_ab_host = ddc::create_mirror_view_and_copy(get_field(collfreq_ab));

        double const me_on_memi(mass(my_ielec) / (mass(my_ielec) + mass(my_iion)));
        ddc::for_each(gridx, [&](IdxX const ix) {
            // test : dlog(T_e - T_i)/dt = -12nu_ei*m_e/(m_e+m_b)
            // should be verified
            double const error = std::fabs(
                    std::log(std::fabs(
                            temperature_host(my_ielec, ix) - temperature_host(my_iion, ix)))
                    - std::log(std::fabs(temperature_init(my_ielec) - temperature_init(my_iion)))
                    + 12 * collfreq_ab_host(my_ielec, ix) * me_on_memi * nbiter * deltat);
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
