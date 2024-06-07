// SPDX-License-Identifier: MIT

#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "chargedensitycalculator.hpp"
#include "constantfluidinitialization.hpp"
#ifdef PERIODIC_RDIMX
#include "femperiodicqnsolver.hpp"
#else
#include "femnonperiodicqnsolver.hpp"
#endif
#include "Lagrange_interpolator.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "nullfluidsolver.hpp"
#include "predcorr.hpp"
#include "predcorr_hybrid.hpp"
#include "qnsolver.hpp"
#include "quadrature.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "spline_interpolator.hpp"
#include "splitrighthandsidesolver.hpp"
#include "splitvlasovsolver.hpp"

/**
 * This test creates one instance of the PredCorr class, and one instance of the 
 * PredCorrFluid class. With the NullFluidSolver class, 
 * this test verifies that the distribution function is the same after applying the 
 * predictor corrector scheme to it, with and without the fluid species.
 */
TEST(GeometryXM, PredCorrHybrid)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_ncells(50);

    CoordVx const vx_min(-8);
    CoordVx const vx_max(8);
    IVectVx const vx_ncells(50);

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
    host_t<FieldSp<int>> fluid_charges(dom_fluidsp);
    ddc::parallel_fill(fluid_charges, 0.);

    host_t<DFieldSp> fluid_masses(dom_fluidsp);
    ddc::parallel_fill(fluid_masses, mass_ion);

    // Create the domain of kinetic species + fluid species
    IDomainSp const dom_allsp(IndexSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains charges of all species
    host_t<FieldSp<int>> charges(dom_allsp);

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

    host_t<FieldSp<int>> init_perturb_mode(dom_kinsp);
    ddc::parallel_fill(init_perturb_mode, 2);
    host_t<DFieldSp> init_perturb_amplitude(dom_kinsp);
    ddc::parallel_fill(init_perturb_amplitude, 0.1);

    SingleModePerturbInitialization const
            init(allfequilibrium,
                 init_perturb_mode.span_cview(),
                 init_perturb_amplitude.span_cview());
    init(allfdistribu);

    // Moments domain initialization
    IVectM const nb_fluid_moments(3);
    IDomainM const meshM(IndexM(0), nb_fluid_moments);
    ddc::init_discrete_space<IDimM>();

    IndexM idensity(0);
    IndexM iflux(1);
    IndexM istress(2);

    // Initialization of fluid species moments
    DFieldSpMX fluid_moments_alloc(IDomainSpMX(dom_fluidsp, meshM, meshX));
    auto fluid_moments = fluid_moments_alloc.span_view();

    host_t<DFieldSpM> moments_init(IDomainSpM(dom_fluidsp, meshM));
    ddc::parallel_fill(moments_init[idensity], 1.);
    ddc::parallel_fill(moments_init[iflux], 0.);
    ddc::parallel_fill(moments_init[istress], 1.);

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

    host_t<DFieldVx> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(meshVx, builder_vx_poisson);

    auto const quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    ChargeDensityCalculator rhs(quadrature_coeffs);
#ifdef PERIODIC_RDIMX
    FFTPoissonSolver<IDomainX, IDomainX, Kokkos::DefaultExecutionSpace> fft_poisson_solver(meshX);
    QNSolver const poisson(fft_poisson_solver, rhs);
#else
    FemNonPeriodicQNSolver const poisson(builder_x_poisson, spline_x_evaluator_poisson, rhs);
#endif

    // construction of predcorr without fluid species
    PredCorr const predcorr(vlasov, poisson);

    // distribution function to be evolved by predcorr without fluid species
    DFieldSpXVx allfdistribu_predcorr_alloc(allfdistribu.domain());
    auto allfdistribu_predcorr = allfdistribu_predcorr_alloc.span_view();
    ddc::parallel_deepcopy(allfdistribu_predcorr, allfdistribu);

    double const time_start(0.);
    int const nb_iter(10);
    double const deltat(0.1);
    predcorr(allfdistribu_predcorr, time_start, deltat, nb_iter);

    // construction of predcorr with fluid species
    NullFluidSolver const fluidsolver(dom_fluidsp);
    PredCorrHybrid const predcorr_hybrid(vlasov, fluidsolver, poisson);
    predcorr_hybrid(allfdistribu, fluid_moments, time_start, deltat, nb_iter);

    auto allfdistribu_host
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), allfdistribu);

    auto allfdistribu_predcorr_host = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), allfdistribu_predcorr);

    /**
     * Since the fluid model uses NullFluidSolver, 
     * the distribution function evolved with PredCorrFluid and PredCorr 
     * should be equal
     */
    double const tolerance(1.e-12);
    ddc::for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        EXPECT_LE(
                std::fabs(allfdistribu_host(ispxvx) - allfdistribu_predcorr_host(ispxvx)),
                tolerance);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
