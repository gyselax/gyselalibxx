// SPDX-License-Identifier: MIT
#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "maxwellianequilibrium.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

static void TestMaxwellianTokamAxi()
{
    /**
     * Tests if the DFieldMemVx maxwellian function defines a maxwellian with the correct moments
     */
    int nb_kinspecies(2);
    IdxRangeSp idx_range_kinsp = IdxRangeSp(IdxSp(0), IdxStepSp(nb_kinspecies));

    CoordR const r_min(0.);
    CoordR const r_max(150.);
    IdxStepR const r_size(4);
    CoordTheta const theta_min(0.);
    CoordTheta const theta_max(2 * M_PI);
    IdxStepTheta const theta_size(4);
    CoordVpar const vpar_min(-7.);
    CoordVpar const vpar_max(7.);
    IdxStepVpar const vpar_size(128);
    CoordMu const mu_min(0.);
    CoordMu const mu_max(12.);
    IdxStepMu const mu_size(256);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesR>(r_min, r_max, r_size);
    ddc::init_discrete_space<BSplinesTheta>(theta_min, theta_max, theta_size);
    ddc::init_discrete_space<BSplinesVpar>(vpar_min, vpar_max, vpar_size);
    ddc::init_discrete_space<BSplinesMu>(mu_min, mu_max, mu_size);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());
    ddc::init_discrete_space<GridVpar>(SplineInterpPointsVpar::get_sampling<GridVpar>());
    ddc::init_discrete_space<GridMu>(SplineInterpPointsMu::get_sampling<GridMu>());

    IdxRangeR gridr(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta gridtheta(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeVpar gridvpar(SplineInterpPointsVpar::get_domain<GridVpar>());
    IdxRangeMu gridmu(SplineInterpPointsMu::get_domain<GridMu>());
    IdxRangeTor2D const idx_range_tor2d(gridr, gridtheta);
    IdxRangeV2D const idx_range_v2d(gridvpar, gridmu);
    IdxRangeV2DTor2D const idx_range_v2dtor2d(gridvpar, gridmu, gridr, gridtheta);
    IdxRangeSpTor2D const idx_range_sptor2d(idx_range_kinsp, gridr, gridtheta);
    IdxRangeSpV2DTor2D const
            idx_range_spv2dtor2d(idx_range_kinsp, gridvpar, gridmu, gridr, gridtheta);

    // Initialization of the density, mean velocity and temperature
    host_t<DFieldMemSpTor2D> density_host_alloc(idx_range_sptor2d);
    host_t<DFieldSpTor2D> density_host = get_field(density_host_alloc);
    host_t<DFieldMemSpTor2D> mean_velocity_host_alloc(idx_range_sptor2d);
    host_t<DFieldSpTor2D> mean_velocity_host = get_field(mean_velocity_host_alloc);
    host_t<DFieldMemSpTor2D> temperature_host_alloc(idx_range_sptor2d);
    host_t<DFieldSpTor2D> temperature_host = get_field(temperature_host_alloc);
    IdxSp idx_sp1 = idx_range_kinsp.front();
    ddc::parallel_fill(Kokkos::DefaultHostExecutionSpace(), density_host[idx_sp1], 1.);
    ddc::parallel_fill(Kokkos::DefaultHostExecutionSpace(), mean_velocity_host[idx_sp1], 0.);
    ddc::parallel_fill(Kokkos::DefaultHostExecutionSpace(), temperature_host[idx_sp1], 1.);
    IdxSp idx_sp2 = idx_range_kinsp.back();
    ddc::parallel_fill(Kokkos::DefaultHostExecutionSpace(), density_host[idx_sp2], 1.2);
    ddc::parallel_fill(Kokkos::DefaultHostExecutionSpace(), mean_velocity_host[idx_sp2], 0.5);
    ddc::parallel_fill(Kokkos::DefaultHostExecutionSpace(), temperature_host[idx_sp2], 1.2);

    // Initialization of the magnetic field
    host_t<DFieldMemTor2D> magnetic_field_host_alloc(idx_range_tor2d);
    host_t<DFieldTor2D> magnetic_field_host = get_field(magnetic_field_host_alloc);
    ddc::parallel_fill(Kokkos::DefaultHostExecutionSpace(), magnetic_field_host, 1.);
    auto magnetic_field_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(magnetic_field_host));
    DFieldTor2D magnetic_field = get_field(magnetic_field_alloc);

    // Initialization of the distribution function as a Maxwellian
    DFieldMemSpV2DTor2D allfequilibrium_alloc(idx_range_spv2dtor2d);
    DFieldSpV2DTor2D allfequilibrium = get_field(allfequilibrium_alloc);
    MaxwellianEquilibrium const init_fequilibrium(
            density_host,
            temperature_host,
            mean_velocity_host,
            magnetic_field_host);
    init_fequilibrium(allfequilibrium);

    // Initialization of the quadrature coefficients for integration in vpar, mu
    DFieldMemV2D quadrature_coeffs_v2d_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(idx_range_v2d));
    Quadrature<IdxRangeV2D, IdxRangeSpV2DTor2D> const integrate_v2d(
            get_const_field(quadrature_coeffs_v2d_alloc));

    // Computation of the density = \int fM Jv dvpar dmu with Jv the normalized velocity Jacobian Jv = 2pi*B
    DFieldMemSpTor2D density_res_alloc(idx_range_sptor2d);
    DFieldSpTor2D density_res = get_field(density_res_alloc);
    integrate_v2d(
            Kokkos::DefaultExecutionSpace(),
            density_res,
            KOKKOS_LAMBDA(IdxSpV2DTor2D const ispv2dtor2d) {
                IdxSp const isp = ddc::select<Species>(ispv2dtor2d);
                IdxTor2D const itor2d = ddc::select<GridR, GridTheta>(ispv2dtor2d);
                double const jacobian_velocity = 2. * M_PI * magnetic_field(itor2d);
                return jacobian_velocity * allfequilibrium(ispv2dtor2d);
            });

    // Computation of the mean velocity Upar
    // --> Computation of n*Upar =  \int fM vpar Jv dvpar dmu
    DFieldMemSpTor2D mean_velocity_res_alloc(idx_range_sptor2d);
    DFieldSpTor2D mean_velocity_res = get_field(mean_velocity_res_alloc);
    integrate_v2d(
            Kokkos::DefaultExecutionSpace(),
            mean_velocity_res,
            KOKKOS_LAMBDA(IdxSpV2DTor2D const ispv2dtor2d) {
                IdxSp const isp = ddc::select<Species>(ispv2dtor2d);
                IdxVpar const ivpar = ddc::select<GridVpar>(ispv2dtor2d);
                IdxTor2D const itor2d = ddc::select<GridR, GridTheta>(ispv2dtor2d);
                double const jacobian_velocity = 2. * M_PI * magnetic_field(itor2d);
                return jacobian_velocity * ddc::coordinate(ivpar) * allfequilibrium(ispv2dtor2d);
            });
    // --> Computation of Upar = n*Upar / n
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_sptor2d,
            KOKKOS_LAMBDA(IdxSpTor2D const isptor2d) {
                mean_velocity_res(isptor2d) = mean_velocity_res(isptor2d) / density_res(isptor2d);
            });

    // Computation of the temperature that is deduced as
    //    T = (Tpar + 2*Tperp)/3
    //  with
    //    Tpar = Ppar/n and Tperp = Pperp/ n (n being the density)
    //  where
    //    . Ppar = \int fM * (vpar-Upar)^2 Jv dvpar dmu
    //    . Perp = \int fM mu B Jv dvpar dmu
    DFieldMemSpTor2D temperature_res_alloc(idx_range_sptor2d);
    DFieldSpTor2D temperature_res = get_field(temperature_res_alloc);
    integrate_v2d(
            Kokkos::DefaultExecutionSpace(),
            temperature_res,
            KOKKOS_LAMBDA(IdxSpV2DTor2D const ispv2dtor2d) {
                IdxSp const isp = ddc::select<Species>(ispv2dtor2d);
                IdxVpar const ivpar = ddc::select<GridVpar>(ispv2dtor2d);
                IdxMu const imu = ddc::select<GridMu>(ispv2dtor2d);
                IdxTor2D const itor2d = ddc::select<GridR, GridTheta>(ispv2dtor2d);
                double const B = magnetic_field(itor2d);
                double const vpar_minus_Upar
                        = ddc::coordinate(ivpar) - mean_velocity_res(isp, itor2d);
                double const mu = ddc::coordinate(imu);
                double const jacobian_velocity = 2. * M_PI * B;
                return jacobian_velocity * (vpar_minus_Upar * vpar_minus_Upar + 2. * mu * B) / 3.
                       * allfequilibrium(ispv2dtor2d);
            });
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_sptor2d,
            KOKKOS_LAMBDA(IdxSpTor2D const isptor2d) {
                temperature_res(isptor2d) = temperature_res(isptor2d) / density_res(isptor2d);
            });

    // Print density, mean velocity and temperature for one (r,theta) value
    IdxRangeR idx_range_tor1_reduce(IdxR(2), IdxStepR(1));
    IdxRangeTheta idx_range_tor2_reduce(IdxTheta(3), IdxStepTheta(1));
    IdxRangeSpTor2D
            idx_range_sptor2d_reduce(idx_range_kinsp, idx_range_tor1_reduce, idx_range_tor2_reduce);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_sptor2d_reduce,
            KOKKOS_LAMBDA(IdxSpTor2D const isptor2d) {
                IdxSp isp = ddc::select<Species>(isptor2d);
                Kokkos::
                        printf("density for sp %d = %e \n",
                               (isp - idx_range_kinsp.front()).value(),
                               density_res(isptor2d));
                Kokkos::
                        printf("mean_velocity for sp %d = %e \n",
                               (isp - idx_range_kinsp.front()).value(),
                               mean_velocity_res(isptor2d));
                Kokkos::
                        printf("temperature for sp %d = %e \n",
                               (isp - idx_range_kinsp.front()).value(),
                               temperature_res(isptor2d));
            });

    // Compute the maximum error for density calculation
    auto density_alloc
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), density_host);
    DFieldSpTor2D density = get_field(density_alloc);
    host_t<DFieldMemSp> max_error_density_alloc(idx_range_kinsp);
    host_t<DFieldSp> max_error_density = get_field(max_error_density_alloc);
    ddc::for_each(idx_range_kinsp, [&](IdxSp const isp) {
        max_error_density(isp) = ddc::parallel_transform_reduce(
                idx_range_tor2d,
                0.0,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IdxTor2D const itor2d) {
                    return Kokkos::fabs(density(isp, itor2d) - density_res(isp, itor2d));
                });
        std::cout << "max_error_density = " << max_error_density(isp) << std::endl;
    });

    // Compute the maximum error for mean velocity calculation
    auto mean_velocity_alloc
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), mean_velocity_host);
    DFieldSpTor2D mean_velocity = get_field(mean_velocity_alloc);
    host_t<DFieldMemSp> max_error_mean_velocity_alloc(idx_range_kinsp);
    host_t<DFieldSp> max_error_mean_velocity = get_field(max_error_mean_velocity_alloc);
    ddc::for_each(idx_range_kinsp, [&](IdxSp const isp) {
        max_error_mean_velocity(isp) = ddc::parallel_transform_reduce(
                idx_range_tor2d,
                0.0,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IdxTor2D const itor2d) {
                    return Kokkos::fabs(
                            mean_velocity(isp, itor2d) - mean_velocity_res(isp, itor2d));
                });
        std::cout << "max_error_mean_velocity = " << max_error_mean_velocity(isp) << std::endl;
    });

    // Compute the maximum error for temperature calculation
    auto temperature_alloc
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), temperature_host);
    DFieldSpTor2D temperature = get_field(temperature_alloc);
    host_t<DFieldMemSp> max_error_temperature_alloc(idx_range_kinsp);
    host_t<DFieldSp> max_error_temperature = get_field(max_error_temperature_alloc);
    ddc::for_each(idx_range_kinsp, [&](IdxSp const isp) {
        max_error_temperature(isp) = ddc::parallel_transform_reduce(
                idx_range_tor2d,
                0.0,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IdxTor2D const itor2d) {
                    return Kokkos::fabs(temperature(isp, itor2d) - temperature_res(isp, itor2d));
                });
        std::cout << "max_error_temperature = " << max_error_temperature(isp) << std::endl;
    });

    ddc::for_each(idx_range_kinsp, [&](IdxSp const isp) {
        EXPECT_LE(max_error_density(isp), 1e-3);
        EXPECT_LE(max_error_mean_velocity(isp), 1e-3);
        EXPECT_LE(max_error_temperature(isp), 1e-3);
    });
}

TEST(MaxwellianTokamAxi, Moments)
{
    TestMaxwellianTokamAxi();
}
