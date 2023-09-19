// SPDX-License-Identifier: MIT

#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <fluid_moments.hpp>
#include <geometry.hpp>
#include <maxwellianequilibrium.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

/**
     * Initializes the distribution function as a Maxwellian with fluid moments depending on space.
     * Uses the FluidMoments methods to recompute these moments.
     */
TEST(Physics, FluidMoments)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_size(10);

    CoordVx const vx_min(-9.);
    CoordVx const vx_max(9.);
    IVectVx const vx_size(400);

    IVectSp const nb_species(2);
    IDomainSp const dom_sp(IndexSp(0), nb_species);
    IndexSp const my_ielec = dom_sp.front();
    IndexSp const my_iion = dom_sp.back();

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
    IDomainSpXVx const mesh(IDomainSp(my_iion, IVectSp(1)), gridx, gridvx);

    FieldSp<int> charges(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    DFieldSp masses(dom_sp);
    ddc::fill(masses, 1);
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
    DFieldSpXVx allfdistribu(mesh);

    // Initialization of the distribution function as a maxwellian with
    // moments depending on space
    DFieldSpX density_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX mean_velocity_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX temperature_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
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
        DFieldVx finit(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                finit.span_view(),
                density_init(ispx),
                temperature_init(ispx),
                mean_velocity_init(ispx));
        ddc::deepcopy(allfdistribu[ispx], finit);
    });

    // density and temperature
    DFieldSpX density_computed(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX mean_velocity_computed(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX temperature_computed(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    FluidMoments moments(Quadrature<IDimVx>(
            trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu))));

    moments(density_computed.span_view(), allfdistribu.span_cview(), FluidMoments::s_density);
    moments(mean_velocity_computed.span_view(),
            allfdistribu.span_cview(),
            density_computed.span_cview(),
            FluidMoments::s_velocity);
    moments(temperature_computed.span_view(),
            allfdistribu.span_cview(),
            density_computed.span_cview(),
            mean_velocity_computed.span_cview(),
            FluidMoments::s_temperature);

    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        EXPECT_LE(std::fabs(density_computed(ispx) - density_init(ispx)), 1e-12);
        EXPECT_LE(std::fabs(mean_velocity_computed(ispx) - mean_velocity_init(ispx)), 1e-12);
        EXPECT_LE(std::fabs(temperature_computed(ispx) - temperature_init(ispx)), 1e-12);
    });
}
