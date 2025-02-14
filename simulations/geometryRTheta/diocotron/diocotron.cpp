// SPDX-License-Identifier: MIT
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_rp.hpp"
#include "bsl_predcorr.hpp"
#include "bsl_predcorr_second_order_explicit.hpp"
#include "bsl_predcorr_second_order_implicit.hpp"
#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "diocotron_initialization_equilibrium.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "l_norm_tools.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "quadrature.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "spline_polar_foot_finder.hpp"
#include "spline_quadrature.hpp"
#include "trapezoid_quadrature.hpp"



namespace {
using PoissonSolver = PolarSplineFEMPoissonLikeSolver<
        GridR,
        GridTheta,
        PolarBSplinesRTheta,
        SplineRThetaEvaluatorNullBound>;
using DiscreteMappingBuilder
        = DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>;
using DiscreteMappingBuilder_host = DiscreteToCartesianBuilder<
        X,
        Y,
        SplineRThetaBuilder_host,
        SplineRThetaEvaluatorConstBound_host>;
// using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
using LogicalToPhysicalMapping = CzarnyToCartesian<R, Theta, X, Y>;

namespace fs = std::filesystem;

} // end namespace



int main(int argc, char** argv)
{
    // SETUP ==========================================================================================
    fs::create_directory("output");

    // Get the parameters of the mesh_rp from the grid_size.yaml. ----------------------------------------
    PC_tree_t conf_gyselalibxx = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;

    start_simulation = std::chrono::system_clock::now();

    // Build the mesh_rp for the space. ------------------------------------------------------------------
    IdxRangeR const mesh_r = init_pseudo_uniform_spline_dependent_idx_range<
            GridR,
            BSplinesR,
            SplineInterpPointsR>(conf_gyselalibxx, "r");
    IdxRangeTheta const mesh_p = init_pseudo_uniform_spline_dependent_idx_range<
            GridTheta,
            BSplinesTheta,
            SplineInterpPointsTheta>(conf_gyselalibxx, "p");
    double const dt(PCpp_double(conf_gyselalibxx, ".Time.delta_t"));
    double const final_T(PCpp_double(conf_gyselalibxx, ".Time.final_T"));

    IdxRangeRTheta const mesh_rp(mesh_r, mesh_p);

    host_t<FieldMemRTheta<CoordRTheta>> coords(mesh_rp);
    ddc::for_each(mesh_rp, [&](IdxRTheta const irp) { coords(irp) = ddc::coordinate(irp); });


    // OPERATORS ======================================================================================
    SplineRThetaBuilder const builder(mesh_rp);
    SplineRThetaBuilder_host const builder_host(mesh_rp);


    // --- Define the mapping. ------------------------------------------------------------------------
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(
            ddc::coordinate(mesh_r.front()));
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(
            ddc::coordinate(mesh_r.back()));

    SplineRThetaEvaluatorConstBound spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());

    SplineRThetaEvaluatorConstBound_host spline_evaluator_extrapol_host(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());

    const LogicalToPhysicalMapping to_physical_mapping(0.3,1.4);
    DiscreteMappingBuilder const discrete_mapping_builder(
            Kokkos::DefaultExecutionSpace(),
            to_physical_mapping,
            builder,
            spline_evaluator_extrapol);
    DiscreteMappingBuilder_host const discrete_mapping_builder_host(
            Kokkos::DefaultHostExecutionSpace(),
            to_physical_mapping,
            builder_host,
            spline_evaluator_extrapol_host);
    DiscreteToCartesian const discrete_mapping = discrete_mapping_builder();
    DiscreteToCartesian const discrete_mapping_host = discrete_mapping_builder_host();


    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping_host);

    IdxRangeBSRTheta const idx_range_bsplinesRTheta = get_spline_idx_range(builder_host);


    // --- Time integration method --------------------------------------------------------------------
#if defined(EULER_METHOD)
    Euler<host_t<FieldMemRTheta<CoordRTheta>>,
          host_t<DVectorFieldMemRTheta<X, Y>>,
          Kokkos::DefaultHostExecutionSpace> const time_stepper(mesh_rp);

#elif defined(CRANK_NICOLSON_METHOD)
    double const epsilon_CN = 1e-8;
    CrankNicolson<
            host_t<FieldMemRTheta<CoordRTheta>>,
            host_t<DVectorFieldMemRTheta<X, Y>>,
            Kokkos::DefaultHostExecutionSpace> const time_stepper(mesh_rp, 20, epsilon_CN);

#elif defined(RK3_METHOD)
    RK3<host_t<FieldMemRTheta<CoordRTheta>>,
        host_t<DVectorFieldMemRTheta<X, Y>>,
        Kokkos::DefaultHostExecutionSpace> const time_stepper(mesh_rp);

#elif defined(RK4_METHOD)
    RK4<host_t<FieldMemRTheta<CoordRTheta>>,
        host_t<DVectorFieldMemRTheta<X, Y>>,
        Kokkos::DefaultHostExecutionSpace> const time_stepper(mesh_rp);

#endif


    // --- Advection operator -------------------------------------------------------------------------
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> p_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);
    SplineRThetaEvaluatorNullBound_host spline_evaluator_host(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    PreallocatableSplineInterpolatorRTheta interpolator(builder_host, spline_evaluator_host);

    SplinePolarFootFinder find_feet(
            time_stepper,
            to_physical_mapping,
            to_physical_mapping,
            builder_host,
            spline_evaluator_extrapol_host);

    BslAdvectionRTheta advection_operator(interpolator, find_feet, to_physical_mapping);



    // --- Poisson solver -----------------------------------------------------------------------------
    // Coefficients alpha and beta of the Poisson equation:
    DFieldMemRTheta coeff_alpha(mesh_rp);
    DFieldMemRTheta coeff_beta(mesh_rp);
    ddc::parallel_fill(coeff_alpha, -1);
    ddc::parallel_fill(coeff_beta, 0);


    Spline2DMem coeff_alpha_spline(idx_range_bsplinesRTheta);
    Spline2DMem coeff_beta_spline(idx_range_bsplinesRTheta);

    builder(get_field(coeff_alpha_spline), get_const_field(coeff_alpha));
    builder(get_field(coeff_beta_spline), get_const_field(coeff_beta));

    PoissonSolver poisson_solver(
            get_const_field(coeff_alpha_spline),
            get_const_field(coeff_beta_spline),
            discrete_mapping,
            spline_evaluator);

    // --- Predictor corrector operator ---------------------------------------------------------------
#if defined(PREDCORR)
    BslPredCorrRTheta predcorr_operator(
            to_physical_mapping,
            advection_operator,
            builder_host,
            spline_evaluator_host,
            poisson_solver);
#elif defined(EXPLICIT_PREDCORR)
    BslExplicitPredCorrRTheta predcorr_operator(
            to_physical_mapping,
            to_physical_mapping,
            advection_operator,
            mesh_rp,
            builder_host,
            spline_evaluator_host,
            poisson_solver,
            spline_evaluator_extrapol_host);
#elif defined(IMPLICIT_PREDCORR)
    BslImplicitPredCorrRTheta predcorr_operator(
            to_physical_mapping,
            to_physical_mapping,
            advection_operator,
            mesh_rp,
            builder_host,
            spline_evaluator_host,
            poisson_solver,
            spline_evaluator_extrapol_host);
#endif

    // ================================================================================================
    // SIMULATION DATA                                                                                 |
    // ================================================================================================
    double const Q(PCpp_double(
            conf_gyselalibxx,
            ".Perturbation.charge_Q")); // no charge carried by the inner conductor r = W1.
    int const l(PCpp_int(conf_gyselalibxx, ".Perturbation.l_mode"));
    double const eps(PCpp_double(conf_gyselalibxx, ".Perturbation.eps"));
    CoordR const R1(PCpp_double(conf_gyselalibxx, ".Perturbation.r_min"));
    CoordR const R2(PCpp_double(conf_gyselalibxx, ".Perturbation.r_max"));
    DiocotronDensitySolution exact_rho(
            ddc::coordinate(mesh_r.front()),
            R1,
            R2,
            ddc::coordinate(mesh_r.back()),
            Q,
            l,
            eps);

    // --- Time parameters ----------------------------------------------------------------------------
    int const iter_nb = final_T * int(1 / dt);

    // --- save simulation data
    ddc::expose_to_pdi("r_size", ddc::discrete_space<BSplinesR>().ncells());
    ddc::expose_to_pdi("p_size", ddc::discrete_space<BSplinesTheta>().ncells());

    expose_mesh_to_pdi("r_coords", mesh_r);
    expose_mesh_to_pdi("p_coords", mesh_p);

    ddc::expose_to_pdi("delta_t", dt);
    ddc::expose_to_pdi("final_T", final_T);
    ddc::expose_to_pdi("time_step_diag", PCpp_int(conf_gyselalibxx, ".Output.time_step_diag"));

    ddc::expose_to_pdi("slope", exact_rho.get_slope());



    // ================================================================================================
    // INITIALISATION                                                                                 |
    // ================================================================================================
    // Cartesian coordinates and jacobian ****************************
    host_t<FieldMemRTheta<CoordX>> coords_x(mesh_rp);
    host_t<FieldMemRTheta<CoordY>> coords_y(mesh_rp);
    host_t<DFieldMemRTheta> jacobian(mesh_rp);
    ddc::for_each(mesh_rp, [&](IdxRTheta const irp) {
        CoordXY coords_xy = to_physical_mapping(ddc::coordinate(irp));
        coords_x(irp) = ddc::select<X>(coords_xy);
        coords_y(irp) = ddc::select<Y>(coords_xy);
        jacobian(irp) = to_physical_mapping.jacobian(ddc::coordinate(irp));
    });



    host_t<DFieldMemRTheta> rho(mesh_rp);
    host_t<DFieldMemRTheta> rho_eq(mesh_rp);

    // Initialize rho and rho equilibrium ****************************
    ddc::for_each(mesh_rp, [&](IdxRTheta const irp) {
        rho(irp) = exact_rho.initialisation(coords(irp));
        rho_eq(irp) = exact_rho.equilibrium(coords(irp));
    });

    // Compute phi equilibrium phi_eq from Poisson solver. ***********
    DFieldMemRTheta phi_eq(mesh_rp);
    host_t<DFieldMemRTheta> phi_eq_host(mesh_rp);
    host_t<Spline2DMem> rho_coef_eq(idx_range_bsplinesRTheta);
    builder_host(get_field(rho_coef_eq), get_const_field(rho_eq));
    PoissonLikeRHSFunction poisson_rhs_eq(get_const_field(rho_coef_eq), spline_evaluator_host);
    poisson_solver(poisson_rhs_eq, get_field(phi_eq));
    ddc::parallel_deepcopy(phi_eq, phi_eq_host);

    // --- Save initial data --------------------------------------------------------------------------
    ddc::PdiEvent("initialization")
            .with("x_coords", coords_x)
            .with("y_coords", coords_y)
            .with("jacobian", jacobian)
            .with("density_eq", rho_eq)
            .with("electrical_potential_eq", phi_eq_host);


    // ================================================================================================
    // SIMULATION                                                                                     |
    // ================================================================================================
    predcorr_operator(get_field(rho), dt, iter_nb);


    end_simulation = std::chrono::system_clock::now();
    display_time_difference("Simulation time: ", start_simulation, end_simulation);


    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
    PC_tree_destroy(&conf_gyselalibxx);

    return EXIT_SUCCESS;
}
