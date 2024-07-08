// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_1d.hpp"
#include "directional_tag.hpp"
#include "euler.hpp"
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "initialization_Kelvin_Helmholtz.hpp"
#include "input.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr_RK2.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_interpolator.hpp"
#include "utils_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_span.hpp"


using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;


int main(int argc, char** argv)
{
    // Create a folder for the output files.
    fs::create_directory("output");

    // Environments variables for profiling
    setenv("KOKKOS_TOOLS_LIBS", KP_KERNEL_TIMER_PATH, false);
    setenv("KOKKOS_TOOLS_TIMER_JSON", "true", false);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    PC_tree_t conf_gyselalibxx = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);


    // CREATING MESH AND SUPPORTS ----------------------------------------------------------------
    IDomainX const interpolation_domain_x = init_spline_dependent_domain<
            IDimX,
            BSplinesX,
            SplineInterpPointsX>(conf_gyselalibxx, "x");

    IDomainY const interpolation_domain_y = init_spline_dependent_domain<
            IDimY,
            BSplinesY,
            SplineInterpPointsY>(conf_gyselalibxx, "y");

    IDomainXY meshXY(interpolation_domain_x, interpolation_domain_y);


    // READING CONFIGURATION ---------------------------------------------------------------------
    // --> Algorithm info
    double const delta_t = PCpp_double(conf_gyselalibxx, ".Algorithm.delta_t");
    double const final_time = PCpp_double(conf_gyselalibxx, ".Algorithm.final_time");
    int const nbiter = int(final_time / delta_t);

    // --> Output info
    int const nbstep_diag = PCpp_int(conf_gyselalibxx, ".Output.nbstep_diag");

    // --> Initial function infos
    double const epsilon = PCpp_double(conf_gyselalibxx, ".PerturbationInfo.perturb_amplitude");
    double const mode_k = 2 * M_PI / ddcHelper::total_interval_length(interpolation_domain_x);


    // DEFINING OPERATORS ------------------------------------------------------------------------
    // Create spline builderss ---
    SplineXBuilder_XY const builder_x(meshXY);
    SplineYBuilder_XY const builder_y(meshXY);

    // Create spline evaluators ---
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
    SplineXEvaluator_XY const spline_x_evaluator(bv_x_min, bv_x_max);

    ddc::PeriodicExtrapolationRule<RDimY> bv_y_min;
    ddc::PeriodicExtrapolationRule<RDimY> bv_y_max;
    SplineYEvaluator_XY const spline_y_evaluator(bv_y_min, bv_y_max);

    // Create spline interpolators ---
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);
    PreallocatableSplineInterpolator const spline_y_interpolator(builder_y, spline_y_evaluator);

    // Create Poisson solver ---
    FFTPoissonSolver<IDomainXY, IDomainXY, Kokkos::DefaultExecutionSpace> const poisson_solver(
            meshXY);

    // Create advection operators ---
    Euler<FieldXY<CoordX>, DFieldXY> euler_x(meshXY);
    BslAdvection1D<
            IDimX,
            IDomainXY,
            IDomainXY,
            SplineXBuilder_XY,
            SplineXEvaluator_XY,
            Euler<FieldXY<CoordX>, DFieldXY>>
            advection_x(spline_x_interpolator, builder_x, spline_x_evaluator, euler_x);

    Euler<FieldXY<CoordY>, DFieldXY> euler_y(meshXY);
    BslAdvection1D<
            IDimY,
            IDomainXY,
            IDomainXY,
            SplineYBuilder_XY,
            SplineYEvaluator_XY,
            Euler<FieldXY<CoordY>, DFieldXY>>
            advection_y(spline_y_interpolator, builder_y, spline_y_evaluator, euler_y);


    // Create an initilizer ---
    KelvinHelmholtzInstabilityInitialization initialize(epsilon, mode_k);


    // Create predcorr operator: predictor-corrector method based on RK2 ---
    PredCorrRK2XY predictor_corrector(poisson_solver, advection_x, advection_y);


    // INITIALISATION ----------------------------------------------------------------------------
    // Initialisation of the distributed function
    DFieldXY allfdistribu_equilibrium_alloc(meshXY);
    DSpanXY allfdistribu_equilibrium = allfdistribu_equilibrium_alloc.span_view();

    DFieldXY allfdistribu_alloc(meshXY);
    DSpanXY allfdistribu = allfdistribu_alloc.span_view();

    initialize(allfdistribu, allfdistribu_equilibrium);


    // Initialisation of the electrostatic potential and electric field (for saving data)
    DFieldXY electrostatic_potential_alloc(meshXY);
    DSpanXY electrostatic_potential = electrostatic_potential_alloc.span_view();

    VectorFieldXY_XY electric_field_alloc(meshXY);
    VectorSpanXY_XY electric_field = electric_field_alloc.span_view();

    poisson_solver(electrostatic_potential, electric_field, allfdistribu);


    // Save the data ---
    /*
        Here, we save 
            * allfdistribu to check the mass conservation; 
            * electric_field to check the energy conservation; 
            * electrostatic_potential, if needed. 
        The data need to be first saved on CPU. 
        The initiliation is saved here but the during the simulation, 
        the data are saved in the predictor-corrector operator. 
    */
    ddc::expose_to_pdi("Nx_spline_cells", PCpp_int(conf_gyselalibxx, ".SplineMesh.x_ncells"));
    ddc::expose_to_pdi("Ny_spline_cells", PCpp_int(conf_gyselalibxx, ".SplineMesh.y_ncells"));
    expose_mesh_to_pdi("MeshX", interpolation_domain_x);
    expose_mesh_to_pdi("MeshY", interpolation_domain_y);

    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("time_step", delta_t);
    ddc::expose_to_pdi("final_time", final_time);

    auto allfdistribu_equilibrium_host = ddc::create_mirror_and_copy(allfdistribu_equilibrium);
    ddc::expose_to_pdi("fdistribu_equilibrium", allfdistribu_equilibrium_host);

    int const iter = 0;
    auto allfdistribu_host = ddc::create_mirror_and_copy(allfdistribu);
    auto electrostatic_potential_host = ddc::create_mirror_and_copy(electrostatic_potential);
    auto electric_field_x_host = ddc::create_mirror_and_copy(ddcHelper::get<RDimX>(electric_field));
    auto electric_field_y_host = ddc::create_mirror_and_copy(ddcHelper::get<RDimY>(electric_field));
    ddc::PdiEvent("iteration")
            .with("iter", iter)
            .and_with("time_saved", iter * delta_t)
            .and_with("fdistribu", allfdistribu_host)
            .and_with("electrostatic_potential", electrostatic_potential_host)
            .and_with("electric_field_x", electric_field_x_host)
            .and_with("electric_field_y", electric_field_y_host);

    // // SIMULATION --------------------------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> const start
            = std::chrono::system_clock::now();

    predictor_corrector(allfdistribu, delta_t, nbiter);

    std::chrono::time_point<std::chrono::system_clock> const end = std::chrono::system_clock::now();

    display_time_difference("Simulation time: ", start, end);

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
    PC_tree_destroy(&conf_gyselalibxx);

    return EXIT_SUCCESS;
}