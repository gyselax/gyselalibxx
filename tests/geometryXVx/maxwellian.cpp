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

static void TestMaxwellian()
{
    /**
     * Tests if the DFieldMemVx maxwellian function defines a maxwellian with the correct moments
     */
    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IdxStepVx const vx_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeVx gridvx(SplineInterpPointsVx::get_domain<GridVx>());

    SplineVxBuilder_1d const builder_vx(gridvx);

    host_t<DFieldMemVx> const quadrature_coeffs_host = trapezoid_quadrature_coefficients(gridvx);
    auto quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(quadrature_coeffs_host));
    Quadrature<IdxRangeVx> integrate_v(get_const_field(quadrature_coeffs));

    double const density = 1.e-5;
    double const mean_velocity = 0.5;
    double const temperature = 0.5;

    DFieldMemVx fdistribu_alloc(gridvx);
    DFieldVx fdistribu = get_field(fdistribu_alloc);
    MaxwellianEquilibrium::compute_maxwellian(fdistribu, density, temperature, mean_velocity);

    // density
    double const density_res = integrate_v(Kokkos::DefaultExecutionSpace(), fdistribu);

    // mean velocity
    DFieldMemVx integrand_alloc(gridvx);
    DFieldVx integrand = get_field(integrand_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvx,
            KOKKOS_LAMBDA(IdxVx const ivx) {
                integrand(ivx) = ddc::coordinate(ivx) * fdistribu(ivx);
            });

    double const mean_velocity_res
            = integrate_v(Kokkos::DefaultExecutionSpace(), get_field(integrand)) / density_res;

    // temperature
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvx,
            KOKKOS_LAMBDA(IdxVx const ivx) {
                double const velocity = ddc::coordinate(ivx) - mean_velocity_res;
                integrand(ivx) = velocity * velocity * fdistribu(ivx);
            });
    double const temperature_res
            = integrate_v(Kokkos::DefaultExecutionSpace(), get_field(integrand)) / density_res;

    EXPECT_LE(std::fabs(density_res - density), 1e-12);
    EXPECT_LE(std::fabs(mean_velocity_res - mean_velocity), 1e-12);
    EXPECT_LE(std::fabs(temperature_res - temperature), 1e-12);
}

TEST(Maxwellian, Moments)
{
    TestMaxwellian();
}