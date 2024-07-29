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
     * Tests if the DFieldVx maxwellian function defines a maxwellian with the correct moments
     */
    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IVectVx const vx_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimVx>(InterpPointsVx::get_sampling());

    IDomainVx gridvx(InterpPointsVx::get_domain());

    SplineVxBuilder_1d const builder_vx(gridvx);

<<<<<<< HEAD
    host_t<DFieldVx> quadrature_coeffs(
            trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(gridvx));
    Quadrature<Kokkos::DefaultHostExecutionSpace, IDimVx> const integrate_v(
            quadrature_coeffs.span_view());
=======
    host_t<DFieldVx> const quadrature_coeffs = trapezoid_quadrature_coefficients(gridvx);
    host_t<Quadrature<IDomainVx>> const integrate_v(quadrature_coeffs);
>>>>>>> origin/main

    double const density = 1.e-5;
    double const mean_velocity = 0.5;
    double const temperature = 0.5;

    host_t<DFieldVx> fdistribu(gridvx);
    MaxwellianEquilibrium::
            compute_maxwellian(fdistribu.span_view(), density, temperature, mean_velocity);

    // density
    double const density_res = integrate_v(fdistribu);

    // mean velocity
    host_t<DFieldVx> integrand(gridvx);
    ddc::for_each(gridvx, [&](IndexVx const ivx) {
        integrand(ivx) = ddc::coordinate(ivx) * fdistribu(ivx);
    });
    double const mean_velocity_res = integrate_v(integrand) / density_res;

    // temperature
    ddc::for_each(gridvx, [&](IndexVx const ivx) {
        double const velocity = ddc::coordinate(ivx) - mean_velocity_res;
        integrand(ivx) = velocity * velocity * fdistribu(ivx);
    });
    double const temperature_res = integrate_v(integrand) / density_res;

    EXPECT_LE(std::fabs(density_res - density), 1e-12);
    EXPECT_LE(std::fabs(mean_velocity_res - mean_velocity), 1e-12);
    EXPECT_LE(std::fabs(temperature_res - temperature), 1e-12);
}
