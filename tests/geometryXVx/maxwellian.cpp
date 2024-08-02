// SPDX-License-Identifier: MIT

#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <geometry.hpp>
#include <maxwellianequilibrium.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

TEST(Maxwellian, Moments)
{
    /**
     * Tests if the DFieldMemVx maxwellian function defines a maxwellian with the correct moments
     */
    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IdxStepVx const vx_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridVx>(InterpPointsVx::get_sampling());

    IdxRangeVx gridvx(InterpPointsVx::get_domain());

    SplineVxBuilder_1d const builder_vx(gridvx);

    host_t<DFieldMemVx> quadrature_coeffs(
            trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(gridvx));
    host_t<Quadrature<IdxRangeVx>> const integrate_v(get_const_field(quadrature_coeffs));

    double const density = 1.e-5;
    double const mean_velocity = 0.5;
    double const temperature = 0.5;

    host_t<DFieldMemVx> fdistribu(gridvx);
    MaxwellianEquilibrium::
            compute_maxwellian(get_field(fdistribu), density, temperature, mean_velocity);

    // density
    double const density_res = integrate_v(fdistribu);

    // mean velocity
    host_t<DFieldMemVx> integrand(gridvx);
    ddc::for_each(gridvx, [&](IdxVx const ivx) {
        integrand(ivx) = ddc::coordinate(ivx) * fdistribu(ivx);
    });
    double const mean_velocity_res = integrate_v(integrand) / density_res;

    // temperature
    ddc::for_each(gridvx, [&](IdxVx const ivx) {
        double const velocity = ddc::coordinate(ivx) - mean_velocity_res;
        integrand(ivx) = velocity * velocity * fdistribu(ivx);
    });
    double const temperature_res = integrate_v(integrand) / density_res;

    EXPECT_LE(std::fabs(density_res - density), 1e-12);
    EXPECT_LE(std::fabs(mean_velocity_res - mean_velocity), 1e-12);
    EXPECT_LE(std::fabs(temperature_res - temperature), 1e-12);
}
