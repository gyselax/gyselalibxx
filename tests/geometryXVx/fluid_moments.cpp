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

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());

    IDomainX interpolation_domain_x(SplineInterpPointsX::get_domain<IDimX>());
    IDomainVx interpolation_domain_vx(SplineInterpPointsVx::get_domain<IDimVx>());

    SplineXBuilder_1d const builder_x(interpolation_domain_x);

    SplineVxBuilder_1d const builder_vx(interpolation_domain_vx);

    IDomainX const gridx = builder_x.interpolation_domain();
    IDomainVx const gridvx = builder_vx.interpolation_domain();
    IDomainSpXVx const mesh(IDomainSp(my_iion, IVectSp(1)), gridx, gridvx);

    host_t<DFieldSp> charges(dom_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldSp> masses(dom_sp);
    ddc::parallel_fill(masses, 1);

    // Initialization of the distribution function as a maxwellian
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));
    host_t<DFieldSpXVx> allfdistribu_host(mesh);
    DFieldSpXVx allfdistribu(mesh);


    // Initialization of the distribution function as a maxwellian with
    // moments depending on space
    host_t<DFieldSpX> density_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    host_t<DFieldSpX> mean_velocity_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    host_t<DFieldSpX> temperature_init(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host), [&](IndexSpX const ispx) {
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

        auto finit_host = ddc::create_mirror_view_and_copy(finit.span_view());
        ddc::parallel_deepcopy(allfdistribu_host[ispx], finit_host);
    });

    // density and temperature
    DFieldSpX density_computed(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    DFieldSpX mean_velocity_computed(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    DFieldSpX temperature_computed(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host));
    host_t<DFieldVx> const quadrature_coeffs = trapezoid_quadrature_coefficients(gridvx);
    Quadrature<IDimVx> integrate(quadrature_coeffs);
    FluidMoments moments(integrate);
    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
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
    auto mean_velocity_computed_host
            = ddc::create_mirror_view_and_copy(mean_velocity_computed.span_view());
    auto temperature_computed_host
            = ddc::create_mirror_view_and_copy(temperature_computed.span_view());
    auto density_computed_host = ddc::create_mirror_view_and_copy(density_computed.span_view());
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu_host), [&](IndexSpX const ispx) {
        EXPECT_LE(std::fabs(density_computed_host(ispx) - density_init(ispx)), 1e-12);
        EXPECT_LE(std::fabs(mean_velocity_computed_host(ispx) - mean_velocity_init(ispx)), 1e-12);
        EXPECT_LE(std::fabs(temperature_computed_host(ispx) - temperature_init(ispx)), 1e-12);
    });
}
