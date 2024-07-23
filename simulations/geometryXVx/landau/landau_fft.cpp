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
#include "chargedensitycalculator.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "maxwellianequilibrium.hpp"
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

    long int iter_start;
    PC_tree_t conf_voicexx;
    parse_executable_arguments(conf_voicexx, iter_start, argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

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

    IDomainSp const dom_kinsp = init_species(conf_voicexx);
    IDomainSpXVx const meshSpXVx(dom_kinsp, meshXVx);
    IDomainSpVx const meshSpVx(dom_kinsp, mesh_vx);

    SplineXBuilder const builder_x(meshXVx);
    SplineVxBuilder const builder_vx(meshXVx);
    SplineVxBuilder_1d const builder_vx_poisson(mesh_vx);

    // Initialization of the distribution function
    DFieldSpVx allfequilibrium(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium
            = MaxwellianEquilibrium::init_from_input(dom_kinsp, conf_voicexx);
    init_fequilibrium(allfequilibrium);

    ddc::expose_to_pdi("iter_start", iter_start);

    DFieldSpXVx allfdistribu(meshSpXVx);
    double time_start(0);
    if (iter_start == 0) {
        SingleModePerturbInitialization const init = SingleModePerturbInitialization::
                init_from_input(allfequilibrium, dom_kinsp, conf_voicexx);
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
    ddc::ConstantExtrapolationRule<RDimX> bv_x_min(ddc::coordinate(mesh_x.front()));
    ddc::ConstantExtrapolationRule<RDimX> bv_x_max(ddc::coordinate(mesh_x.back()));
#endif

    // Creating operators
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    ddc::ConstantExtrapolationRule<RDimVx> bv_v_min(ddc::coordinate(mesh_vx.front()));
    ddc::ConstantExtrapolationRule<RDimVx> bv_v_max(ddc::coordinate(mesh_vx.back()));

    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);
    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    BslAdvectionSpatial<GeometryXVx, IDimX> const advection_x(spline_x_interpolator);
    BslAdvectionVelocity<GeometryXVx, IDimVx> const advection_vx(spline_vx_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    DFieldVx const quadrature_coeffs(neumann_spline_quadrature_coefficients<
                                     Kokkos::DefaultExecutionSpace>(mesh_vx, builder_vx_poisson));

    ChargeDensityCalculator rhs(quadrature_coeffs.span_view());
    FFTPoissonSolver<IDomainX, IDomainX, Kokkos::DefaultExecutionSpace> fft_poisson_solver(mesh_x);
    QNSolver const poisson(fft_poisson_solver, rhs);

    PredCorr const predcorr(vlasov, poisson);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
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
