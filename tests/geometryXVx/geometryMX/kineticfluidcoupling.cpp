// SPDX-License-Identifier: MIT

#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

#include "Lagrange_interpolator.hpp"
#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "charge_exchange.hpp"
#include "chargedensitycalculator.hpp"
#include "constantfluidinitialization.hpp"
#include "diffusiveneutralsolver.hpp"
#include "fem_1d_poisson_solver.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "ionization.hpp"
#include "irighthandside.hpp"
#include "kinetic_fluid_coupling_source.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "predcorr.hpp"
#include "predcorr_hybrid.hpp"
#include "qnsolver.hpp"
#include "quadrature.hpp"
#include "recombination.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "spline_interpolator.hpp"
#include "splitrighthandsidesolver.hpp"
#include "splitvlasovsolver.hpp"

/**
 * This test initializes the fluid species with constant reaction rates for ionization and recombination.
 * Then, using the analytical solution for this scenario where T is cte. we compare it to the solver.
 */
static void TestKineticFluidCoupling()
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_ncells(10);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IVectVx const vx_ncells(20);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_ncells);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());

    IDomainX meshX(SplineInterpPointsX::get_domain<IDimX>());
    IDomainVx meshVx(SplineInterpPointsVx::get_domain<IDimVx>());
    IDomainXVx meshXVx(meshX, meshVx);

    SplineXBuilder const builder_x(meshXVx);
#ifndef PERIODIC_RDIMX
    SplineXBuilder_1d const builder_x_poisson(meshX);
#endif
    SplineVxBuilder const builder_vx(meshXVx);
    SplineVxBuilder_1d const builder_vx_poisson(meshVx);

    // Kinetic species domain initialization
    IVectSp const nb_kinspecies(2);
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IndexSp const iion = dom_kinsp.front();
    IndexSp const ielec = dom_kinsp.back();

    host_t<FieldSp<int>> kinetic_charges(dom_kinsp);
    kinetic_charges(ielec) = -1;
    kinetic_charges(iion) = 1;

    host_t<DFieldSp> kinetic_masses(dom_kinsp);
    double const mass_ion(400), mass_elec(1);
    kinetic_masses(ielec) = mass_elec;
    kinetic_masses(iion) = mass_ion;

    // Fluid species domain initialization
    IVectSp const nb_fluidspecies(1);
    IDomainSp const dom_fluidsp(IndexSp(dom_kinsp.back() + 1), nb_fluidspecies);

    // Fluid charges
    host_t<DFieldSp> fluid_charges(dom_fluidsp);
    ddc::parallel_fill(fluid_charges, 0.);

    host_t<DFieldSp> fluid_masses(dom_fluidsp);
    ddc::parallel_fill(fluid_masses, mass_ion);

    // Create the domain of kinetic species + fluid species
    IDomainSp const dom_allsp(IndexSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains charges of all species
    host_t<DFieldSp> charges(dom_allsp);

    // fill the Field with charges of kinetic species
    for (IndexSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }

    // fill the Field with charges of fluid species
    for (IndexSp isp : dom_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }

    // Create a Field that contains masses of kinetic and fluid species
    host_t<DFieldSp> masses(dom_allsp);

    // fill the Field with masses of kinetic species
    for (IndexSp isp : dom_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }

    // fill the Field with masses of fluid species
    for (IndexSp isp : dom_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

    // Initialization of kinetic species distribution function
    DFieldSpXVx allfdistribu_alloc(IDomainSpXVx(dom_kinsp, meshX, meshVx));
    auto allfdistribu = allfdistribu_alloc.span_view();

    host_t<DFieldSp> kinsp_density_eq(dom_kinsp);
    host_t<DFieldSp> kinsp_velocity_eq(dom_kinsp);
    host_t<DFieldSp> kinsp_temperature_eq(dom_kinsp);

    ddc::parallel_fill(kinsp_density_eq, 1.);
    ddc::parallel_fill(kinsp_velocity_eq, 0.);
    ddc::parallel_fill(kinsp_temperature_eq, 1.);

    DFieldSpVx allfequilibrium_alloc(IDomainSpVx(dom_kinsp, meshVx));
    auto allfequilibrium = allfequilibrium_alloc.span_view();
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(kinsp_density_eq),
            std::move(kinsp_temperature_eq),
            std::move(kinsp_velocity_eq));
    init_fequilibrium(allfequilibrium);

    host_t<IFieldSp> init_perturb_mode(dom_kinsp);
    host_t<DFieldSp> init_perturb_amplitude(dom_kinsp);
    ddc::parallel_fill(init_perturb_mode, 1);
    ddc::parallel_fill(init_perturb_amplitude, 0.0);

    SingleModePerturbInitialization const
            init(allfequilibrium, std::move(init_perturb_mode), std::move(init_perturb_amplitude));
    init(allfdistribu);

    // Moments domain initialization
    IVectM const nb_fluid_moments(1);
    IDomainM const meshM(IndexM(0), nb_fluid_moments);
    ddc::init_discrete_space<IDimM>();

    // Initialization of fluid species moments
    DFieldSpMX fluid_moments_alloc(IDomainSpMX(dom_fluidsp, meshM, meshX));
    auto fluid_moments = fluid_moments_alloc.span_view();

    host_t<DFieldSpM> moments_init(IDomainSpM(dom_fluidsp, meshM));
    ddc::parallel_fill(moments_init, 0.);
    ConstantFluidInitialization fluid_init(moments_init);
    fluid_init(fluid_moments);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
#else
    ddc::ConstantExtrapolationRule<RDimX> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<RDimX> bv_x_max(x_max);
#endif

    // Creating operators
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
#ifndef PERIODIC_RDIMX
    SplineXEvaluator_1d const spline_x_evaluator_poisson(bv_x_min, bv_x_max);
#endif
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    IVectVx static constexpr gwvx {0};
    LagrangeInterpolator<IDimVx, BCond::DIRICHLET, BCond::DIRICHLET, IDimX, IDimVx> const
            lagrange_vx_non_preallocatable_interpolator(3, gwvx);
    PreallocatableLagrangeInterpolator<
            IDimVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            IDimX,
            IDimVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);

    BslAdvectionSpatial<GeometryXVx, IDimX> const advection_x(spline_x_interpolator);
    BslAdvectionVelocity<GeometryXVx, IDimVx> const advection_vx(lagrange_vx_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    DFieldVx const quadrature_coeffs(neumann_spline_quadrature_coefficients<
                                     Kokkos::DefaultExecutionSpace>(meshVx, builder_vx_poisson));

    ChargeDensityCalculator rhs(quadrature_coeffs.span_cview());
#ifdef PERIODIC_RDIMX
    FFTPoissonSolver<IDomainX, IDomainX, Kokkos::DefaultExecutionSpace> poisson_solver(meshX);
#else
    FEM1DPoissonSolver const poisson_solver(builder_x_poisson, spline_x_evaluator_poisson);
#endif
    QNSolver const poisson(poisson_solver, rhs);

    double const normalization_coeff(0.01);
    double const k_0(1.e-3);

    ChargeExchangeRate charge_exchange(k_0);
    IonizationRate ionization(k_0);
    RecombinationRate recombination(k_0);

    SplineXBuilder_1d const spline_x_builder_neutrals(meshX);
    SplineXEvaluator_1d const spline_x_evaluator_neutrals(bv_x_min, bv_x_max);

    DFieldVx const quadrature_coeffs_neutrals(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(meshVx));

    DiffusiveNeutralSolver const fluidsolver(
            charge_exchange,
            ionization,
            recombination,
            normalization_coeff,
            spline_x_builder_neutrals,
            spline_x_evaluator_neutrals,
            quadrature_coeffs_neutrals.span_cview());

    // kinetic fluid coupling term
    KineticFluidCouplingSource const kineticfluidcoupling(
            1.,
            0.,
            0.,
            ionization,
            recombination,
            normalization_coeff,
            quadrature_coeffs);

    double const time_start(0.);
    int const nb_iter(20);
    double const deltat(0.1);

    PredCorrHybrid const predcorr_hybrid(vlasov, fluidsolver, poisson, kineticfluidcoupling);
    predcorr_hybrid(allfdistribu, fluid_moments, time_start, deltat, nb_iter);

    auto allfdistribu_host
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), allfdistribu);
    auto fluid_moments_host
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), fluid_moments);

    // analytical solution
    // we know the rate values for the initial conditions, we assume T=cte. so rates independent of time.
    double const ionization_rate = 1.130359390036803 * k_0;
    double const recombination_rate = 7.638123065868132e-06 * k_0;

    double const N = 1.;
    double const alpha = -(N * ionization_rate) / normalization_coeff;
    double const beta = -(ionization_rate + recombination_rate) / (2 * normalization_coeff);
    double const C = (recombination_rate + ionization_rate)
                     / (2 * (0.0 * ionization_rate - 1.0 * recombination_rate));
    double const X_1
            = N * (recombination_rate - ionization_rate) / (recombination_rate + ionization_rate);

    DFieldSpMX X_alloc(IDomainSpMX(dom_fluidsp, meshM, meshX));
    auto X = X_alloc.span_view();
    DFieldSpMX analytical_nN_alloc(IDomainSpMX(dom_fluidsp, meshM, meshX));
    auto analytical_nN = analytical_nN_alloc.span_view();
    double const t_diag = nb_iter * deltat;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            X.domain(),
            KOKKOS_LAMBDA(IndexSpMX const ispmx) {
                X(ispmx) = Kokkos::exp(alpha * t_diag)
                                   * Kokkos::
                                           pow((beta / alpha) * (Kokkos::exp(alpha * t_diag) - 1)
                                                       + C,
                                               -1)
                           + X_1;
                analytical_nN(ispmx) = (X(ispmx) + N) / 2.;
            });

    auto analytical_nN_host = ddc::create_mirror_view_and_copy(analytical_nN);

    ddc::for_each(fluid_moments_host.domain(), [&](IndexSpMX const ispmx) {
        EXPECT_NEAR(analytical_nN_host(ispmx), fluid_moments_host(ispmx), 1.5e-8);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}

TEST(GeometryMX, KineticFluidCoupling)
{
    TestKineticFluidCoupling();
}
