// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <collisions_inter.hpp>
#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "collisions_intra.hpp"
#include "constantfluidinitialization.hpp"
#include "constantrate.hpp"
#include "diffusiveneutralsolver.hpp"
#ifdef PERIODIC_RDIMX
#include "femperiodicqnsolver.hpp"
#else
#include "femnonperiodicqnsolver.hpp"
#endif
#include "Lagrange_interpolator.hpp"
#include "chargedensitycalculator.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "irighthandside.hpp"
#include "kinetic_source.hpp"
#include "krook_source_adaptive.hpp"
#include "krook_source_constant.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "neutrals.yml.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "pdi_out_neutrals.yml.hpp"
#include "predcorr_hybrid.hpp"
#include "qnsolver.hpp"
#include "restartinitialization.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "species_init.hpp"
#include "spline_interpolator.hpp"
#include "splitrighthandsidesolver.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    // Environments variables for profiling
    setenv("KOKKOS_TOOLS_LIBS", KP_KERNEL_TIMER_PATH, false);
    setenv("KOKKOS_TOOLS_TIMER_JSON", "true", false);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    long int iter_start;
    PC_tree_t conf_voicexx;
    parse_executable_arguments(conf_voicexx, iter_start, argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    // Reading config
    // --> Mesh info
    IDomainX const mesh_x = init_spline_dependent_domain<
            IDimX,
            BSplinesX,
            SplineInterpPointsX>(conf_voicexx, "x");
    IDomainVx const mesh_vx = init_spline_dependent_domain<
            IDimVx,
            BSplinesVx,
            SplineInterpPointsVx>(conf_voicexx, "vx");
    IDomainXVx const meshXVx(mesh_x, mesh_vx);

    SplineXBuilder const builder_x(meshXVx);
#ifndef PERIODIC_RDIMX
    SplineXBuilder_1d const builder_x_poisson(mesh_x);
#endif
    SplineVxBuilder const builder_vx(meshXVx);
    SplineVxBuilder_1d const builder_vx_poisson(mesh_vx);

    IDomainSp dom_kinsp;
    IDomainSp dom_fluidsp;
    init_species_withfluid(dom_kinsp, dom_fluidsp, conf_voicexx);

    host_t<DFieldSp> kinetic_density_eq(dom_kinsp);
    host_t<DFieldSp> kinetic_temperature_eq(dom_kinsp);
    host_t<DFieldSp> kinetic_mean_velocity_eq(dom_kinsp);

    host_t<DFieldSp> init_perturb_amplitude(dom_kinsp);
    host_t<FieldSp<int>> init_perturb_mode(dom_kinsp);

    for (IndexSp const isp : dom_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(conf_voicexx, ".SpeciesInfo[%d]", isp.uid());

        kinetic_density_eq(isp) = PCpp_double(conf_isp, ".density_eq");
        kinetic_temperature_eq(isp) = PCpp_double(conf_isp, ".temperature_eq");
        kinetic_mean_velocity_eq(isp) = PCpp_double(conf_isp, ".mean_velocity_eq");
        init_perturb_amplitude(isp) = PCpp_double(conf_isp, ".perturb_amplitude");
        init_perturb_mode(isp) = static_cast<int>(PCpp_int(conf_isp, ".perturb_mode"));
    }

    // Initialization of kinetic species distribution function
    IDomainSpVx const meshSpVx(dom_kinsp, mesh_vx);
    DFieldSpVx allfequilibrium(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(kinetic_density_eq),
            std::move(kinetic_temperature_eq),
            std::move(kinetic_mean_velocity_eq));
    init_fequilibrium(allfequilibrium);

    ddc::expose_to_pdi("iter_start", iter_start);

    IDomainSpXVx const meshSpXVx(dom_kinsp, meshXVx);
    DFieldSpXVx allfdistribu(meshSpXVx);
    double time_start(0);
    if (iter_start == 0) {
        SingleModePerturbInitialization const
                init(allfequilibrium,
                     init_perturb_mode.span_cview(),
                     init_perturb_amplitude.span_cview());
        init(allfdistribu);
    } else {
        RestartInitialization const restart(iter_start, time_start);
        restart(allfdistribu);
    }
    auto allfequilibrium_host = ddc::create_mirror_view_and_copy(allfequilibrium.span_view());

    // Moments domain initialization
    IVectM const nb_fluid_moments(1);
    IDomainM const meshM(IndexM(0), nb_fluid_moments);
    ddc::init_discrete_space<IDimM>();

    // Neutral species initialization
    DFieldSpMX neutrals_alloc(IDomainSpMX(dom_fluidsp, meshM, mesh_x));
    auto neutrals = neutrals_alloc.span_view();

    host_t<DFieldSpM> moments_init(IDomainSpM(dom_fluidsp, meshM));
    ddc::parallel_fill(moments_init, 1.);
    ConstantFluidInitialization fluid_init(moments_init);
    fluid_init(neutrals);

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
#else
    ddc::ConstantExtrapolationRule<RDimX> bv_x_min(ddc::coordinate(mesh_x.front()));
    ddc::ConstantExtrapolationRule<RDimX> bv_x_max(ddc::coordinate(mesh_x.back()));
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

    // list of rhs operators
    std::vector<std::reference_wrapper<IRightHandSide const>> rhs_operators;
    std::vector<KrookSourceConstant> krook_source_constant_vector;
    std::vector<KrookSourceAdaptive> krook_source_adaptive_vector;
    // Krook operators initialization
    int const nb_rhsKrook(PCpp_len(conf_voicexx, ".Krook"));
    for (int ik = 0; ik < nb_rhsKrook; ++ik) {
        // --> Krook info
        PC_tree_t const conf_krook = PCpp_get(conf_voicexx, ".Krook[%d]", ik);

        static std::map<std::string, RhsType>
                str2rhstype {{"source", RhsType::Source}, {"sink", RhsType::Sink}};
        RhsType type = str2rhstype[PCpp_string(conf_krook, ".type")];
        std::string const krook_name = PCpp_string(conf_krook, ".name");
        if (krook_name == "constant") {
            krook_source_constant_vector.emplace_back(
                    mesh_x,
                    mesh_vx,
                    type,
                    PCpp_double(conf_krook, ".extent"),
                    PCpp_double(conf_krook, ".stiffness"),
                    PCpp_double(conf_krook, ".amplitude"),
                    PCpp_double(conf_krook, ".density"),
                    PCpp_double(conf_krook, ".temperature"));
            rhs_operators.emplace_back(krook_source_constant_vector.back());

        } else if (krook_name == "adaptive") {
            krook_source_adaptive_vector.emplace_back(
                    mesh_x,
                    mesh_vx,
                    type,
                    PCpp_double(conf_krook, ".extent"),
                    PCpp_double(conf_krook, ".stiffness"),
                    PCpp_double(conf_krook, ".amplitude"),
                    PCpp_double(conf_krook, ".density"),
                    PCpp_double(conf_krook, ".temperature"));
            rhs_operators.emplace_back(krook_source_adaptive_vector.back());
        } else {
            throw std::invalid_argument(
                    "Invalid krook name, allowed values are: 'constant', or 'adaptive'.");
        }
    }

    // Kinetic source
    KineticSource const rhs_kinetic_source(
            mesh_x,
            mesh_vx,
            PCpp_double(conf_voicexx, ".KineticSource.extent"),
            PCpp_double(conf_voicexx, ".KineticSource.stiffness"),
            PCpp_double(conf_voicexx, ".KineticSource.amplitude"),
            PCpp_double(conf_voicexx, ".KineticSource.density"),
            PCpp_double(conf_voicexx, ".KineticSource.energy"),
            PCpp_double(conf_voicexx, ".KineticSource.temperature"));
    rhs_operators.emplace_back(rhs_kinetic_source);


    CollisionsIntra const
            collisions_intra(meshSpXVx, PCpp_double(conf_voicexx, ".CollisionsInfo.nustar0"));
    rhs_operators.emplace_back(collisions_intra);

    std::optional<CollisionsInter> collisions_inter;
    if (PCpp_bool(conf_voicexx, ".CollisionsInfo.enable_inter")) {
        collisions_inter.emplace(meshSpXVx, PCpp_double(conf_voicexx, ".CollisionsInfo.nustar0"));
        rhs_operators.emplace_back(*collisions_inter);
    }
    SplitVlasovSolver const vlasov(advection_x, advection_vx);
    SplitRightHandSideSolver const boltzmann(vlasov, rhs_operators);

    host_t<DFieldVx> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(mesh_vx, builder_vx_poisson);

    auto const quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    ChargeDensityCalculator rhs(quadrature_coeffs);
#ifdef PERIODIC_RDIMX
    FFTPoissonSolver<IDomainX, IDomainX, Kokkos::DefaultExecutionSpace> fft_poisson_solver(mesh_x);
    QNSolver const poisson(fft_poisson_solver, rhs);
#else
    FemNonPeriodicQNSolver const poisson(builder_x_poisson, spline_x_evaluator_poisson, rhs);
#endif

    double const charge_exchange_val(1.);
    double const ionization_val(1.e-6);
    double const recombination_val(1.e-7);
    ConstantRate charge_exchange(charge_exchange_val);
    ConstantRate ionization(ionization_val);
    ConstantRate recombination(recombination_val);

    double const neutrals_temperature(1.);
    double const normalization_coeff(0.01);

    SplineXBuilder_1d const spline_x_builder_neutrals(mesh_x);
    SplineXEvaluator_1d const spline_x_evaluator_neutrals(bv_x_min, bv_x_max);

    host_t<DFieldVx> const quadrature_coeffs_neutrals_host(
            trapezoid_quadrature_coefficients(mesh_vx));
    auto const quadrature_coeffs_neutrals = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_neutrals_host.span_view());

    DiffusiveNeutralSolver const neutralsolver(
            charge_exchange,
            ionization,
            recombination,
            neutrals_temperature,
            normalization_coeff,
            spline_x_builder_neutrals,
            spline_x_evaluator_neutrals,
            quadrature_coeffs_neutrals.span_cview());

    PredCorrHybrid const predcorr(boltzmann, neutralsolver, poisson);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
    expose_mesh_to_pdi("MeshX", mesh_x);
    expose_mesh_to_pdi("MeshVx", mesh_vx);
    ddc::expose_to_pdi("Lx", ddcHelper::total_interval_length(mesh_x));
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", dom_kinsp.size());
    ddc::expose_to_pdi("fdistribu_charges", ddc::discrete_space<IDimSp>().charges()[dom_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<IDimSp>().masses()[dom_kinsp]);
    ddc::expose_to_pdi("neutrals_masses", ddc::discrete_space<IDimSp>().masses()[dom_fluidsp]);
    ddc::expose_to_pdi("charge_exchange", charge_exchange_val);
    ddc::expose_to_pdi("ionization", ionization_val);
    ddc::expose_to_pdi("recombination", recombination_val);
    ddc::expose_to_pdi("normalization_coeff_neutrals", normalization_coeff);
    ddc::expose_to_pdi("neutrals_temperature", neutrals_temperature);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu, neutrals, time_start, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
