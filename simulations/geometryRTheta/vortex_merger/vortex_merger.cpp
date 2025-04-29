// SPDX-License-Identifier: MIT
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_polar.hpp"
#include "bsl_predcorr.hpp"
#include "bsl_predcorr_second_order_explicit.hpp"
#include "bsl_predcorr_second_order_implicit.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "l_norm_tools.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "quadrature.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_interpolator_2d.hpp"
#include "spline_polar_foot_finder.hpp"
#include "spline_quadrature.hpp"
#include "trapezoid_quadrature.hpp"
#include "vortex_merger_equilibrium.hpp"
#include "vortex_merger_initialisation.hpp"


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
using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;

} // end namespace

namespace fs = std::filesystem;


int main(int argc, char** argv)
{
    // SETUP ==========================================================================================
    fs::create_directory("output");

    // Get the parameters of the grid from the grid_size.yaml. ----------------------------------------
    PC_tree_t conf_gyselalibxx = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;

    start_simulation = std::chrono::system_clock::now();

    // Build the grid for the space. ------------------------------------------------------------------
    int const Nr(PCpp_int(conf_gyselalibxx, ".SplineMesh.r_ncells"));
    int const Nt(PCpp_int(conf_gyselalibxx, ".SplineMesh.theta_ncells"));
    double const dt(PCpp_double(conf_gyselalibxx, ".Time.delta_t"));
    double const final_T(PCpp_double(conf_gyselalibxx, ".Time.final_T"));

    IdxRangeR const mesh_r = init_pseudo_uniform_spline_dependent_idx_range<
            GridR,
            BSplinesR,
            SplineInterpPointsR>(conf_gyselalibxx, "r");
    IdxRangeTheta const mesh_theta = init_pseudo_uniform_spline_dependent_idx_range<
            GridTheta,
            BSplinesTheta,
            SplineInterpPointsTheta>(conf_gyselalibxx, "theta");
    IdxRangeRTheta const grid(mesh_r, mesh_theta);


    // OPERATORS ======================================================================================
    SplineRThetaBuilder const builder(grid);
    SplineRThetaBuilder_host const builder_host(grid);

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

    const LogicalToPhysicalMapping to_physical_mapping;
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
    Euler<FieldMemRTheta<CoordRTheta>,
          DVectorFieldMemRTheta<X, Y>,
          Kokkos::DefaultExecutionSpace> const time_stepper(grid);


    // --- Advection operator -------------------------------------------------------------------------
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);
    SplineRThetaEvaluatorNullBound_host spline_evaluator_host(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

    PreallocatableSplineInterpolator2D interpolator(builder, spline_evaluator, grid);

    SplinePolarFootFinder find_feet(
            time_stepper,
            to_physical_mapping,
            to_physical_mapping,
            builder,
            spline_evaluator_extrapol);

    BslAdvectionPolar advection_operator(interpolator, find_feet, to_physical_mapping);



    // --- Poisson solver -----------------------------------------------------------------------------
    // Coefficients alpha and beta of the Poisson equation:
    DFieldMemRTheta coeff_alpha(grid);
    DFieldMemRTheta coeff_beta(grid);

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
    BslImplicitPredCorrRTheta predcorr_operator(
            to_physical_mapping,
            to_physical_mapping,
            advection_operator,
            grid,
            builder_host,
            poisson_solver,
            spline_evaluator_extrapol_host);



    // ================================================================================================
    // SIMULATION DATA                                                                                 |
    // ================================================================================================
    double const eps(PCpp_double(conf_gyselalibxx, ".Perturbation.eps"));
    const double sigma(PCpp_double(conf_gyselalibxx, ".Perturbation.sigma"));
    const double x_star_1(PCpp_double(conf_gyselalibxx, ".Perturbation.x_star_1"));
    const double y_star_1(PCpp_double(conf_gyselalibxx, ".Perturbation.y_star_1"));
    const double x_star_2(PCpp_double(conf_gyselalibxx, ".Perturbation.x_star_2"));
    const double y_star_2(PCpp_double(conf_gyselalibxx, ".Perturbation.y_star_2"));

    // --- Time parameters ----------------------------------------------------------------------------
    int const iter_nb = final_T * int(1 / dt);

    // --- save simulation data
    ddc::expose_to_pdi("r_size", Nr);
    ddc::expose_to_pdi("theta_size", Nt);

    host_t<FieldMemR<CoordR>> coords_r(ddc::select<GridR>(grid));
    host_t<FieldMemTheta<CoordTheta>> coords_p(ddc::select<GridTheta>(grid));
    ddc::for_each(ddc::select<GridR>(grid), [&](IdxR const ir) {
        coords_r(ir) = ddc::coordinate(ir);
    });
    ddc::for_each(ddc::select<GridTheta>(grid), [&](IdxTheta const itheta) {
        coords_p(itheta) = ddc::coordinate(itheta);
    });

    ddc::expose_to_pdi("r_coords", coords_r);
    ddc::expose_to_pdi("theta_coords", coords_p);

    ddc::expose_to_pdi("delta_t", dt);
    ddc::expose_to_pdi("final_T", final_T);
    ddc::expose_to_pdi("time_step_diag", PCpp_int(conf_gyselalibxx, ".Output.time_step_diag"));


    // ================================================================================================
    // INITIALISATION                                                                                 |
    // ================================================================================================
    // Cartesian coordinates and jacobian ****************************
    host_t<FieldMemRTheta<CoordX>> coords_x(grid);
    host_t<FieldMemRTheta<CoordY>> coords_y(grid);
    host_t<DFieldMemRTheta> jacobian(grid);
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        CoordXY coords_xy = to_physical_mapping(ddc::coordinate(irtheta));
        coords_x(irtheta) = ddc::select<X>(coords_xy);
        coords_y(irtheta) = ddc::select<Y>(coords_xy);
        jacobian(irtheta) = to_physical_mapping.jacobian(ddc::coordinate(irtheta));
    });

    double const tau(1e-10);
    double const phi_max(1.);


    VortexMergerEquilibria equilibrium(
            to_physical_mapping,
            grid,
            builder_host,
            spline_evaluator_host,
            poisson_solver);
    std::function<double(double const)> const function = [&](double const x) { return x * x; };
    host_t<DFieldMemRTheta> rho_eq(grid);
    equilibrium.set_equilibrium(get_field(rho_eq), function, phi_max, tau);


    VortexMergerDensitySolution solution(to_physical_mapping);
    host_t<DFieldMemRTheta> rho(grid);
    solution.set_initialisation(
            get_field(rho),
            get_const_field(rho_eq),
            eps,
            sigma,
            x_star_1,
            y_star_1,
            x_star_2,
            y_star_2);


    // Compute phi equilibrium phi_eq from Poisson solver. ***********
    DFieldMemRTheta phi_eq(grid);
    host_t<DFieldMemRTheta> phi_eq_host(grid);
    host_t<Spline2DMem> rho_coef_eq(idx_range_bsplinesRTheta);
    builder_host(get_field(rho_coef_eq), get_const_field(rho_eq));
    PoissonLikeRHSFunction poisson_rhs_eq(get_const_field(rho_coef_eq), spline_evaluator_host);
    poisson_solver(poisson_rhs_eq, get_field(phi_eq));
    ddc::parallel_deepcopy(phi_eq, phi_eq_host);


    // --- Save initial data --------------------------------------------------------------------------
    ddc::PdiEvent("initialisation")
            .with("x_coords", coords_x)
            .with("y_coords", coords_y)
            .with("jacobian", jacobian)
            .with("density_eq", rho_eq)
            .with("electrical_potential_eq", phi_eq);



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
