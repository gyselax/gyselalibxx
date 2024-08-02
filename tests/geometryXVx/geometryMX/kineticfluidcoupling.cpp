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
    IdxStepX const x_ncells(10);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IdxStepVx const vx_ncells(20);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_ncells);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX meshX(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx meshVx(SplineInterpPointsVx::get_domain<GridVx>());
    IdxRangeXVx meshXVx(meshX, meshVx);

    SplineXBuilder const builder_x(meshXVx);
#ifndef PERIODIC_RDIMX
    SplineXBuilder_1d const builder_x_poisson(meshX);
#endif
    SplineVxBuilder const builder_vx(meshXVx);
    SplineVxBuilder_1d const builder_vx_poisson(meshVx);

    // Kinetic species index range initialization
    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const dom_kinsp(IdxSp(0), nb_kinspecies);

    IdxSp const iion = dom_kinsp.front();
    IdxSp const ielec = dom_kinsp.back();

    host_t<FieldMemSp<int>> kinetic_charges(dom_kinsp);
    kinetic_charges(ielec) = -1;
    kinetic_charges(iion) = 1;

    host_t<DFieldMemSp> kinetic_masses(dom_kinsp);
    double const mass_ion(400), mass_elec(1);
    kinetic_masses(ielec) = mass_elec;
    kinetic_masses(iion) = mass_ion;

    // Fluid species index range initialization
    IdxStepSp const nb_fluidspecies(1);
    IdxRangeSp const dom_fluidsp(IdxSp(dom_kinsp.back() + 1), nb_fluidspecies);

    // Fluid charges
    host_t<DFieldMemSp> fluid_charges(dom_fluidsp);
    ddc::parallel_fill(fluid_charges, 0.);

    host_t<DFieldMemSp> fluid_masses(dom_fluidsp);
    ddc::parallel_fill(fluid_masses, mass_ion);

    // Create the index range of kinetic species + fluid species
    IdxRangeSp const dom_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains charges of all species
    host_t<DFieldMemSp> charges(dom_allsp);

    // fill the Field with charges of kinetic species
    for (IdxSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }

    // fill the Field with charges of fluid species
    for (IdxSp isp : dom_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }

    // Create a Field that contains masses of kinetic and fluid species
    host_t<DFieldMemSp> masses(dom_allsp);

    // fill the Field with masses of kinetic species
    for (IdxSp isp : dom_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }

    // fill the Field with masses of fluid species
    for (IdxSp isp : dom_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    // Initialization of kinetic species distribution function
    DFieldMemSpXVx allfdistribu_alloc(IdxRangeSpXVx(dom_kinsp, meshX, meshVx));
    auto allfdistribu = get_field(allfdistribu_alloc);

    host_t<DFieldMemSp> kinsp_density_eq(dom_kinsp);
    host_t<DFieldMemSp> kinsp_velocity_eq(dom_kinsp);
    host_t<DFieldMemSp> kinsp_temperature_eq(dom_kinsp);

    ddc::parallel_fill(kinsp_density_eq, 1.);
    ddc::parallel_fill(kinsp_velocity_eq, 0.);
    ddc::parallel_fill(kinsp_temperature_eq, 1.);

    DFieldMemSpVx allfequilibrium_alloc(IdxRangeSpVx(dom_kinsp, meshVx));
    auto allfequilibrium = get_field(allfequilibrium_alloc);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(kinsp_density_eq),
            std::move(kinsp_temperature_eq),
            std::move(kinsp_velocity_eq));
    init_fequilibrium(allfequilibrium);

    host_t<IFieldMemSp> init_perturb_mode(dom_kinsp);
    host_t<DFieldMemSp> init_perturb_amplitude(dom_kinsp);
    ddc::parallel_fill(init_perturb_mode, 1);
    ddc::parallel_fill(init_perturb_amplitude, 0.0);

    SingleModePerturbInitialization const
            init(allfequilibrium, std::move(init_perturb_mode), std::move(init_perturb_amplitude));
    init(allfdistribu);

    // Moments index range initialization
    IdxStepMom const nb_fluid_moments(1);
    IdxRangeMom const meshM(IdxMom(0), nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();

    // Initialization of fluid species moments
    DFieldMemSpMomX fluid_moments_alloc(IdxRangeSpMomX(dom_fluidsp, meshM, meshX));
    auto fluid_moments = get_field(fluid_moments_alloc);

    host_t<DFieldMemSpMom> moments_init(IdxRangeSpMom(dom_fluidsp, meshM));
    ddc::parallel_fill(moments_init, 0.);
    ConstantFluidInitialization fluid_init(moments_init);
    fluid_init(fluid_moments);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
#else
    ddc::ConstantExtrapolationRule<X> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<X> bv_x_max(x_max);
#endif

    // Creating operators
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
#ifndef PERIODIC_RDIMX
    SplineXEvaluator_1d const spline_x_evaluator_poisson(bv_x_min, bv_x_max);
#endif
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    IdxStepVx static constexpr gwvx {0};
    LagrangeInterpolator<GridVx, BCond::DIRICHLET, BCond::DIRICHLET, GridX, GridVx> const
            lagrange_vx_non_preallocatable_interpolator(3, gwvx);
    PreallocatableLagrangeInterpolator<
            GridVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            GridX,
            GridVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);

    BslAdvectionSpatial<GeometryXVx, GridX> const advection_x(spline_x_interpolator);
    BslAdvectionVelocity<GeometryXVx, GridVx> const advection_vx(lagrange_vx_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    host_t<DFieldMemVx> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(meshVx, builder_vx_poisson);

    auto const quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(quadrature_coeffs_host));
    ChargeDensityCalculator rhs(quadrature_coeffs);
#ifdef PERIODIC_RDIMX
    FFTPoissonSolver<IdxRangeX, IdxRangeX, Kokkos::DefaultExecutionSpace> poisson_solver(meshX);
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

    host_t<DFieldMemVx> const quadrature_coeffs_neutrals_host(
            trapezoid_quadrature_coefficients(meshVx));
    auto const quadrature_coeffs_neutrals = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(quadrature_coeffs_neutrals_host));

    DiffusiveNeutralSolver const fluidsolver(
            charge_exchange,
            ionization,
            recombination,
            normalization_coeff,
            spline_x_builder_neutrals,
            spline_x_evaluator_neutrals,
            get_const_field(quadrature_coeffs_neutrals));

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

    DFieldMemSpMomX X_alloc(IdxRangeSpMomX(dom_fluidsp, meshM, meshX));
    auto X = get_field(X_alloc);
    DFieldMemSpMomX analytical_nN_alloc(IdxRangeSpMomX(dom_fluidsp, meshM, meshX));
    auto analytical_nN = get_field(analytical_nN_alloc);
    double const t_diag = nb_iter * deltat;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(X),
            KOKKOS_LAMBDA(IdxSpMomX const ispmx) {
                X(ispmx) = Kokkos::exp(alpha * t_diag)
                                   * Kokkos::
                                           pow((beta / alpha) * (Kokkos::exp(alpha * t_diag) - 1)
                                                       + C,
                                               -1)
                           + X_1;
                analytical_nN(ispmx) = (X(ispmx) + N) / 2.;
            });

    auto analytical_nN_host = ddc::create_mirror_view_and_copy(analytical_nN);

    ddc::for_each(get_idx_range(fluid_moments_host), [&](IdxSpMomX const ispmx) {
        EXPECT_NEAR(analytical_nN_host(ispmx), fluid_moments_host(ispmx), 1.5e-8);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}

TEST(GeometryMX, KineticFluidCoupling)
{
    TestKineticFluidCoupling();
}
