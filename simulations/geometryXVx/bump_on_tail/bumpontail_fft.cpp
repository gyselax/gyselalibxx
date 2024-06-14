// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "bumpontailequilibrium.hpp"
#include "chargedensitycalculator.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "neumann_spline_quadrature.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "qnsolver.hpp"
#include "restartinitialization.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "species_init.hpp"
#include "spline_interpolator.hpp"
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
    CoordX const x_min(PCpp_double(conf_voicexx, ".SplineMesh.x_min"));
    CoordX const x_max(PCpp_double(conf_voicexx, ".SplineMesh.x_max"));
    IVectX const x_ncells(PCpp_int(conf_voicexx, ".SplineMesh.x_ncells"));
    CoordVx const vx_min(PCpp_double(conf_voicexx, ".SplineMesh.vx_min"));
    CoordVx const vx_max(PCpp_double(conf_voicexx, ".SplineMesh.vx_max"));
    IVectVx const vx_ncells(PCpp_int(conf_voicexx, ".SplineMesh.vx_ncells"));

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_ncells);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    IDomainX mesh_x(SplineInterpPointsX::get_domain<IDimX>());
    IDomainVx mesh_vx(SplineInterpPointsVx::get_domain<IDimVx>());
    IDomainXVx meshXVx(mesh_x, mesh_vx);

    IDomainSp const dom_kinsp = init_species(conf_voicexx);

    IDomainSpXVx const meshSpXVx(dom_kinsp, mesh_x, mesh_vx);
    IDomainSpVx const meshSpVx(dom_kinsp, mesh_vx);

    SplineXBuilder const builder_x(meshXVx);
    SplineVxBuilder const builder_vx(meshXVx);
    SplineVxBuilder_1d const builder_vx_poisson(mesh_vx);

    host_t<DFieldSp> epsilon_bot(dom_kinsp);
    host_t<DFieldSp> temperature_bot(dom_kinsp);
    host_t<DFieldSp> mean_velocity_bot(dom_kinsp);
    host_t<DFieldSp> init_perturb_amplitude(dom_kinsp);
    host_t<FieldSp<int>> init_perturb_mode(dom_kinsp);

    for (IndexSp const isp : dom_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(conf_voicexx, ".SpeciesInfo[%d]", isp.uid());

        epsilon_bot(isp) = PCpp_double(conf_isp, ".epsilon_bot");
        temperature_bot(isp) = PCpp_double(conf_isp, ".temperature_bot");
        mean_velocity_bot(isp) = PCpp_double(conf_isp, ".mean_velocity_bot");
        init_perturb_amplitude(isp) = PCpp_double(conf_isp, ".perturb_amplitude");
        init_perturb_mode(isp) = static_cast<int>(PCpp_int(conf_isp, ".perturb_mode"));
    }

    // Initialization of the distribution function
    DFieldSpVx allfequilibrium(meshSpVx);
    BumpontailEquilibrium const init_fequilibrium(
            std::move(epsilon_bot),
            std::move(temperature_bot),
            std::move(mean_velocity_bot));
    init_fequilibrium(allfequilibrium);

    ddc::expose_to_pdi("iter_start", iter_start);

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
    ddc::ConstantExtrapolationRule<RDimX> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<RDimX> bv_x_max(x_max);
#endif

    // Creating operators
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    ddc::ConstantExtrapolationRule<RDimVx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<RDimVx> bv_v_max(vx_max);

    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);
    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    BslAdvectionSpatial<GeometryXVx, IDimX> const advection_x(spline_x_interpolator);
    BslAdvectionVelocity<GeometryXVx, IDimVx> const advection_vx(spline_vx_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    host_t<DFieldVx> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(mesh_vx, builder_vx_poisson);
    auto const quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    FFTPoissonSolver<IDomainX, IDomainX, Kokkos::DefaultExecutionSpace> fft_poisson_solver(mesh_x);
    ChargeDensityCalculator rhs(quadrature_coeffs);
    QNSolver const poisson(fft_poisson_solver, rhs);

    PredCorr const predcorr(vlasov, poisson);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", x_ncells.value());
    ddc::expose_to_pdi("Nvx_spline_cells", vx_ncells.value());
    expose_mesh_to_pdi("MeshX", mesh_x);
    expose_mesh_to_pdi("MeshVx", mesh_vx);
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", dom_kinsp.size());
    ddc::expose_to_pdi("fdistribu_charges", ddc::discrete_space<IDimSp>().charges()[dom_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<IDimSp>().masses()[dom_kinsp]);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu, time_start, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
