// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
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

    PC_tree_t conf_voicexx = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    // Reading config
    // --> Mesh info
    IdxRangeX const mesh_x = init_spline_dependent_idx_range<
            GridX,
            BSplinesX,
            SplineInterpPointsX>(conf_voicexx, "x");
    IdxRangeY const mesh_y = init_spline_dependent_idx_range<
            GridY,
            BSplinesY,
            SplineInterpPointsY>(conf_voicexx, "y");
    IdxRangeVx const mesh_vx = init_spline_dependent_idx_range<
            GridVx,
            BSplinesVx,
            SplineInterpPointsVx>(conf_voicexx, "vx");
    IdxRangeVy const mesh_vy = init_spline_dependent_idx_range<
            GridVy,
            BSplinesVy,
            SplineInterpPointsVy>(conf_voicexx, "vy");
    IdxRangeXY const mesh_xy(mesh_x, mesh_y);
    IdxRangeVxVy mesh_vxvy(mesh_vx, mesh_vy);
    IdxRangeXYVxVy const meshXYVxVy(mesh_x, mesh_y, mesh_vx, mesh_vy);

    IdxRangeSp const dom_kinsp = init_species(conf_voicexx);

    IdxRangeSpVxVy const meshSpVxVy(dom_kinsp, mesh_vx, mesh_vy);
    IdxRangeSpXYVxVy const meshSpXYVxVy(dom_kinsp, meshXYVxVy);

    SplineXBuilder const builder_x(meshXYVxVy);
    SplineYBuilder const builder_y(meshXYVxVy);
    SplineVxBuilder const builder_vx(meshXYVxVy);
    SplineVyBuilder const builder_vy(meshXYVxVy);
    SplineVxBuilder_1d const builder_vx_1d(mesh_vx);
    SplineVyBuilder_1d const builder_vy_1d(mesh_vy);

    // Initialization of the distribution function
    DFieldMemSpVxVy allfequilibrium(meshSpVxVy);
    MaxwellianEquilibrium const init_fequilibrium
            = MaxwellianEquilibrium::init_from_input(dom_kinsp, conf_voicexx);
    init_fequilibrium(allfequilibrium);
    DFieldMemSpXYVxVy allfdistribu(meshSpXYVxVy);
    SingleModePerturbInitialization const init = SingleModePerturbInitialization::
            init_from_input(allfequilibrium, dom_kinsp, conf_voicexx);
    init(allfdistribu);
    auto allfequilibrium_host = ddc::create_mirror_view_and_copy(get_field(allfequilibrium));

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    // Create spline evaluator
    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    ddc::PeriodicExtrapolationRule<Y> bv_y_min;
    ddc::PeriodicExtrapolationRule<Y> bv_y_max;
    SplineYEvaluator const spline_y_evaluator(bv_y_min, bv_y_max);

    PreallocatableSplineInterpolator const spline_y_interpolator(builder_y, spline_y_evaluator);

    ddc::ConstantExtrapolationRule<Vx> bv_vx_min(ddc::coordinate(mesh_vx.front()));
    ddc::ConstantExtrapolationRule<Vx> bv_vx_max(ddc::coordinate(mesh_vx.back()));
    SplineVxEvaluator const spline_vx_evaluator(bv_vx_min, bv_vx_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    ddc::ConstantExtrapolationRule<Vy> bv_vy_min(ddc::coordinate(mesh_vy.front()));
    ddc::ConstantExtrapolationRule<Vy> bv_vy_max(ddc::coordinate(mesh_vy.back()));
    SplineVyEvaluator const spline_vy_evaluator(bv_vy_min, bv_vy_max);

    PreallocatableSplineInterpolator const spline_vy_interpolator(builder_vy, spline_vy_evaluator);

    // Create advection operator
    BslAdvectionSpatial<GeometryXYVxVy, GridX> const advection_x(spline_x_interpolator);
    BslAdvectionSpatial<GeometryXYVxVy, GridY> const advection_y(spline_y_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, GridVx> const advection_vx(spline_vx_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, GridVy> const advection_vy(spline_vy_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_y, advection_vx, advection_vy);

    DFieldMemVxVy const quadrature_coeffs(
            neumann_spline_quadrature_coefficients<
                    Kokkos::DefaultExecutionSpace>(mesh_vxvy, builder_vx_1d, builder_vy_1d));

    FFTPoissonSolver<IDomainXY, IDomainXY, Kokkos::DefaultExecutionSpace> fft_poisson_solver(
            mesh_xy);
    ChargeDensityCalculator const rhs(get_const(quadrature_coeffs));
    QNSolver const poisson(fft_poisson_solver, rhs);

    // Create predcorr operator
    PredCorr const predcorr(vlasov, poisson);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    ddc::expose_to_pdi("Ny_spline_cells", ddc::discrete_space<BSplinesY>().ncells());
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
    ddc::expose_to_pdi("Nvy_spline_cells", ddc::discrete_space<BSplinesVy>().ncells());
    expose_mesh_to_pdi("MeshX", mesh_x);
    expose_mesh_to_pdi("MeshY", mesh_y);
    expose_mesh_to_pdi("MeshVx", mesh_vx);
    expose_mesh_to_pdi("MeshVy", mesh_vy);
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", dom_kinsp.size());
    ddc::expose_to_pdi("fdistribu_charges", ddc::discrete_space<Species>().charges()[dom_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<Species>().masses()[dom_kinsp]);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(get_field(allfdistribu), deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
