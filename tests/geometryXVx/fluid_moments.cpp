// SPDX-License-Identifier: MIT
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "fluid_moments.hpp"
#include "geometry.hpp"
#include "maxwellianequilibrium.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

/**
     * Initialises the distribution function as a Maxwellian with fluid moments depending on space.
     * Uses the FluidMoments methods to recompute these moments.
     */
TEST(Physics, FluidMoments)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(10);

    CoordVx const vx_min(-9.);
    CoordVx const vx_max(9.);
    IdxStepVx const vx_size(400);

    IdxStepSp const nb_species(2);
    IdxRangeSp const idx_range_sp(IdxSp(0), nb_species);
    IdxSp const my_ielec = idx_range_sp.front();
    IdxSp const my_iion = idx_range_sp.back();

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx gridvx(SplineInterpPointsVx::get_domain<GridVx>());

    SplineXBuilder_1d const builder_x(gridx);
    SplineVxBuilder_1d const builder_vx(gridvx);

    IdxRangeSpXVx const mesh(IdxRangeSp(my_iion, IdxStepSp(1)), gridx, gridvx);

    host_t<DFieldMemSp> charges(idx_range_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    host_t<DFieldMemSp> masses(idx_range_sp);
    ddc::parallel_fill(masses, 1);

    // Initialisation of the distribution function as a maxwellian
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
    host_t<DFieldMemSpXVx> allfdistribu_host(mesh);
    DFieldMemSpXVx allfdistribu(mesh);


    // Initialisation of the distribution function as a maxwellian with
    // moments depending on space
    host_t<DFieldMemSpX> density_init(get_idx_range<Species, GridX>(allfdistribu_host));
    host_t<DFieldMemSpX> mean_velocity_init(get_idx_range<Species, GridX>(allfdistribu_host));
    host_t<DFieldMemSpX> temperature_init(get_idx_range<Species, GridX>(allfdistribu_host));
    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu_host), [&](IdxSpX const ispx) {
        double const density = 1.;
        double const density_ampl = 0.1;
        double const mean_velocity = 0.;
        double const mean_velocity_ampl = 0.2;
        double const temperature = 1;
        double const temperature_ampl = 0.3;

        double const coordx = ddc::coordinate(ddc::select<GridX>(ispx));
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
        DFieldMemVx finit(gridvx);
        MaxwellianEquilibrium::compute_maxwellian(
                get_field(finit),
                density_init(ispx),
                temperature_init(ispx),
                mean_velocity_init(ispx));

        auto finit_host = ddc::create_mirror_view_and_copy(get_field(finit));
        ddc::parallel_deepcopy(allfdistribu_host[ispx], finit_host);
    });

    // density and temperature
    DFieldMemSpX density_computed(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldMemSpX mean_velocity_computed(get_idx_range<Species, GridX>(allfdistribu_host));
    DFieldMemSpX temperature_computed(get_idx_range<Species, GridX>(allfdistribu_host));

    DFieldMemVx const quadrature_coeffs
            = trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridvx);
    Quadrature<IdxRangeVx, IdxRangeSpXVx> integrate(get_const_field(quadrature_coeffs));

    FluidMoments moments(integrate);
    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    moments(get_field(density_computed), get_const_field(allfdistribu), FluidMoments::s_density);
    moments(get_field(mean_velocity_computed),
            get_const_field(allfdistribu),
            get_const_field(density_computed),
            FluidMoments::s_velocity);
    moments(get_field(temperature_computed),
            get_const_field(allfdistribu),
            get_const_field(density_computed),
            get_const_field(mean_velocity_computed),
            FluidMoments::s_temperature);
    auto mean_velocity_computed_host
            = ddc::create_mirror_view_and_copy(get_field(mean_velocity_computed));
    auto temperature_computed_host
            = ddc::create_mirror_view_and_copy(get_field(temperature_computed));
    auto density_computed_host = ddc::create_mirror_view_and_copy(get_field(density_computed));
    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu_host), [&](IdxSpX const ispx) {
        EXPECT_LE(std::fabs(density_computed_host(ispx) - density_init(ispx)), 1e-12);
        EXPECT_LE(std::fabs(mean_velocity_computed_host(ispx) - mean_velocity_init(ispx)), 1e-12);
        EXPECT_LE(std::fabs(temperature_computed_host(ispx) - temperature_init(ispx)), 1e-12);
    });
}
