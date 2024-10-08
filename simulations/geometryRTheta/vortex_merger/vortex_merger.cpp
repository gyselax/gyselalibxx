// SPDX-License-Identifier: MIT
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_builder.hpp>
#include <sll/mapping/discrete_to_cartesian.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "advection_domain.hpp"
#include "bsl_advection_rp.hpp"
#include "bsl_predcorr.hpp"
#include "bsl_predcorr_second_order_explicit.hpp"
#include "bsl_predcorr_second_order_implicit.hpp"
#include "compute_norms.hpp"
#include "crank_nicolson.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "poisson_like_rhs_function.hpp"
#include "polarpoissonlikesolver.hpp"
#include "quadrature.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "spline_quadrature.hpp"
#include "trapezoid_quadrature.hpp"
#include "utils_tools.hpp"
#include "vortex_merger_equilibrium.hpp"
#include "vortex_merger_initialization.hpp"


namespace {
using PoissonSolver = PolarSplineFEMPoissonLikeSolver;
using DiscreteMappingBuilder
        = DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>;
using CircularMapping = CircularToCartesian<X, Y, R, Theta>;

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
    int const Nr(PCpp_int(conf_gyselalibxx, ".Mesh.r_size"));
    int const Nt(PCpp_int(conf_gyselalibxx, ".Mesh.p_size"));
    double const dt(PCpp_double(conf_gyselalibxx, ".Time.delta_t"));
    double const final_T(PCpp_double(conf_gyselalibxx, ".Time.final_T"));

    IdxRangeR const mesh_r = init_pseudo_uniform_spline_dependent_idx_range<
            GridR,
            BSplinesR,
            SplineInterpPointsR>(conf_gyselalibxx, "r");
    IdxRangeTheta const mesh_p = init_pseudo_uniform_spline_dependent_idx_range<
            GridTheta,
            BSplinesTheta,
            SplineInterpPointsTheta>(conf_gyselalibxx, "p");
    IdxRangeRTheta const grid(mesh_r, mesh_p);

    host_t<FieldMemRTheta<CoordRTheta>> coords(grid);
    ddc::for_each(grid, [&](IdxRTheta const irp) { coords(irp) = ddc::coordinate(irp); });


    // OPERATORS ======================================================================================
    SplineRThetaBuilder const builder(grid);

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

    const CircularMapping mapping;
    DiscreteMappingBuilder const discrete_mapping_builder(
            Kokkos::DefaultHostExecutionSpace(),
            mapping,
            builder,
            spline_evaluator_extrapol);
    DiscreteToCartesian const discrete_mapping = discrete_mapping_builder();

    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping);

    IdxRangeBSRTheta const idx_range_bsplinesRTheta = get_spline_idx_range(builder);


    // --- Time integration method --------------------------------------------------------------------
    Euler<host_t<FieldMemRTheta<CoordRTheta>>, host_t<DVectorFieldMemRTheta<X, Y>>> const
            time_stepper(grid);


    // --- Advection operator -------------------------------------------------------------------------
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> p_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    PreallocatableSplineInterpolatorRTheta interpolator(builder, spline_evaluator);

    AdvectionPhysicalDomain advection_idx_range(mapping);

    SplineFootFinder
            find_feet(time_stepper, advection_idx_range, builder, spline_evaluator_extrapol);

    BslAdvectionRTheta advection_operator(interpolator, find_feet, mapping);



    // --- Poisson solver -----------------------------------------------------------------------------
    // Coefficients alpha and beta of the Poisson equation:
    host_t<DFieldMemRTheta> coeff_alpha(grid);
    host_t<DFieldMemRTheta> coeff_beta(grid);

    ddc::for_each(grid, [&](IdxRTheta const irp) {
        coeff_alpha(irp) = -1.0;
        coeff_beta(irp) = 0.0;
    });

    host_t<Spline2DMem> coeff_alpha_spline(idx_range_bsplinesRTheta);
    host_t<Spline2DMem> coeff_beta_spline(idx_range_bsplinesRTheta);

    builder(get_field(coeff_alpha_spline), get_const_field(coeff_alpha));
    builder(get_field(coeff_beta_spline), get_const_field(coeff_beta));

    PoissonSolver poisson_solver(coeff_alpha_spline, coeff_beta_spline, discrete_mapping);

    // --- Predictor corrector operator ---------------------------------------------------------------
    BslImplicitPredCorrRTheta predcorr_operator(
            advection_idx_range,
            mapping,
            advection_operator,
            grid,
            builder,
            spline_evaluator,
            poisson_solver,
            spline_evaluator_extrapol);



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
    ddc::expose_to_pdi("p_size", Nt);

    host_t<FieldMemR<CoordR>> coords_r(ddc::select<GridR>(grid));
    host_t<FieldMemTheta<CoordTheta>> coords_p(ddc::select<GridTheta>(grid));
    ddc::for_each(ddc::select<GridR>(grid), [&](IdxR const ir) {
        coords_r(ir) = ddc::coordinate(ir);
    });
    ddc::for_each(ddc::select<GridTheta>(grid), [&](IdxTheta const ip) {
        coords_p(ip) = ddc::coordinate(ip);
    });

    ddc::expose_to_pdi("r_coords", coords_r);
    ddc::expose_to_pdi("p_coords", coords_p);

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
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        CoordXY coords_xy = mapping(ddc::coordinate(irp));
        coords_x(irp) = ddc::select<X>(coords_xy);
        coords_y(irp) = ddc::select<Y>(coords_xy);
        jacobian(irp) = mapping.jacobian(ddc::coordinate(irp));
    });

    double const tau(1e-10);
    double const phi_max(1.);


    VortexMergerEquilibria equilibrium(mapping, grid, builder, spline_evaluator, poisson_solver);
    std::function<double(double const)> const function = [&](double const x) { return x * x; };
    host_t<DFieldMemRTheta> rho_eq(grid);
    equilibrium.set_equilibrium(rho_eq, function, phi_max, tau);


    VortexMergerDensitySolution solution(mapping);
    host_t<DFieldMemRTheta> rho(grid);
    solution.set_initialisation(rho, rho_eq, eps, sigma, x_star_1, y_star_1, x_star_2, y_star_2);


    // Compute phi equilibrium phi_eq from Poisson solver. ***********
    host_t<DFieldMemRTheta> phi_eq(grid);
    host_t<Spline2DMem> rho_coef_eq(idx_range_bsplinesRTheta);
    builder(get_field(rho_coef_eq), get_const_field(rho_eq));
    PoissonLikeRHSFunction poisson_rhs_eq(get_field(rho_coef_eq), spline_evaluator);
    poisson_solver(poisson_rhs_eq, get_const_field(coords), get_field(phi_eq));


    // --- Save initial data --------------------------------------------------------------------------
    ddc::PdiEvent("initialization")
            .with("x_coords", coords_x)
            .and_with("y_coords", coords_y)
            .and_with("jacobian", jacobian)
            .and_with("density_eq", rho_eq)
            .and_with("electrical_potential_eq", phi_eq);



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
