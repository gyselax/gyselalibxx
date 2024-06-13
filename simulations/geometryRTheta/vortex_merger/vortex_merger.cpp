// SPDX-License-Identifier: MIT

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>

#include <paraconf.h>
#include <pdi.h>
#include <utils_tools.hpp>

#include "advection_domain.hpp"
#include "bsl_advection_rp.hpp"
#include "bsl_predcorr.hpp"
#include "bsl_predcorr_second_order_explicit.hpp"
#include "bsl_predcorr_second_order_implicit.hpp"
#include "compute_norms.hpp"
#include "crank_nicolson.hpp"
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
#include "vortex_merger_equilibrium.hpp"
#include "vortex_merger_initialization.hpp"


namespace {
using PoissonSolver = PolarSplineFEMPoissonLikeSolver;
using DiscreteMapping
        = DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder, SplineRPEvaluatorConstBound>;
using CircularMapping = CircularToCartesian<RDimX, RDimY, RDimR, RDimP>;

} // end namespace

namespace fs = std::filesystem;


int main(int argc, char** argv)
{
    // SETUP ==========================================================================================
    fs::create_directory("output");

    // Get the parameters of the grid from the grid_size.yaml. ----------------------------------------
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);

    PC_tree_t conf_gyselalibxx;
    if (argc == 2) {
        conf_gyselalibxx = PC_parse_path(fs::path(argv[1]).c_str());
    } else if (argc == 3) {
        if (argv[1] == std::string_view("--dump-config")) {
            std::fstream file(argv[2], std::fstream::out);
            file << params_yaml;
            return EXIT_SUCCESS;
        }
    } else {
        std::cerr << "usage: " << argv[0] << " [--dump-config] <config_file.yml>" << std::endl;
        return EXIT_FAILURE;
    }
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;

    start_simulation = std::chrono::system_clock::now();

    // Build the grid for the space. ------------------------------------------------------------------
    int const Nr(PCpp_int(conf_gyselalibxx, ".Mesh.r_size"));
    int const Nt(PCpp_int(conf_gyselalibxx, ".Mesh.p_size"));
    double const dt(PCpp_double(conf_gyselalibxx, ".Time.delta_t"));
    double const final_T(PCpp_double(conf_gyselalibxx, ".Time.final_T"));

    IDomainR const mesh_r = init_pseudo_uniform_spline_dependent_domain<
            IDimR,
            BSplinesR,
            SplineInterpPointsR>(conf_gyselalibxx, "r");
    IDomainP const mesh_p = init_pseudo_uniform_spline_dependent_domain<
            IDimP,
            BSplinesP,
            SplineInterpPointsP>(conf_gyselalibxx, "p");
    IDomainRP const grid(mesh_r, mesh_p);

    FieldRP<CoordRP> coords(grid);
    ddc::for_each(grid, [&](IndexRP const irp) { coords(irp) = ddc::coordinate(irp); });


    // OPERATORS ======================================================================================
    SplineRPBuilder const builder(grid);

    // --- Define the mapping. ------------------------------------------------------------------------
    ddc::ConstantExtrapolationRule<RDimR, RDimP> boundary_condition_r_left(
            ddc::coordinate(mesh_r.front()));
    ddc::ConstantExtrapolationRule<RDimR, RDimP> boundary_condition_r_right(
            ddc::coordinate(mesh_r.back()));

    SplineRPEvaluatorConstBound spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<RDimP>(),
            ddc::PeriodicExtrapolationRule<RDimP>());

    const CircularMapping mapping;
    DiscreteMapping const discrete_mapping
            = DiscreteMapping::analytical_to_discrete(mapping, builder, spline_evaluator_extrapol);

    ddc::init_discrete_space<PolarBSplinesRP>(discrete_mapping);

    BSDomainRP const dom_bsplinesRP = builder.spline_domain();


    // --- Time integration method --------------------------------------------------------------------
    Euler<FieldRP<CoordRP>, VectorDFieldRP<RDimX, RDimY>> const time_stepper(grid);


    // --- Advection operator -------------------------------------------------------------------------
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<RDimP> p_extrapolation_rule;
    SplineRPEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    PreallocatableSplineInterpolatorRP interpolator(builder, spline_evaluator);

    AdvectionPhysicalDomain advection_domain(mapping);

    SplineFootFinder find_feet(time_stepper, advection_domain, builder, spline_evaluator_extrapol);

    BslAdvectionRP advection_operator(interpolator, find_feet, mapping);



    // --- Poisson solver -----------------------------------------------------------------------------
    // Coefficients alpha and beta of the Poisson equation:
    DFieldRP coeff_alpha(grid);
    DFieldRP coeff_beta(grid);

    ddc::for_each(grid, [&](IndexRP const irp) {
        coeff_alpha(irp) = -1.0;
        coeff_beta(irp) = 0.0;
    });

    Spline2D coeff_alpha_spline(dom_bsplinesRP);
    Spline2D coeff_beta_spline(dom_bsplinesRP);

    builder(coeff_alpha_spline.span_view(), coeff_alpha.span_cview());
    builder(coeff_beta_spline.span_view(), coeff_beta.span_cview());

    PoissonSolver poisson_solver(coeff_alpha_spline, coeff_beta_spline, discrete_mapping);

    // --- Predictor corrector operator ---------------------------------------------------------------
    BslImplicitPredCorrRP predcorr_operator(
            advection_domain,
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

    FieldR<CoordR> coords_r(ddc::select<IDimR>(grid));
    FieldP<CoordP> coords_p(ddc::select<IDimP>(grid));
    ddc::for_each(ddc::select<IDimR>(grid), [&](IndexR const ir) {
        coords_r(ir) = ddc::coordinate(ir);
    });
    ddc::for_each(ddc::select<IDimP>(grid), [&](IndexP const ip) {
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
    FieldRP<CoordX> coords_x(grid);
    FieldRP<CoordY> coords_y(grid);
    DFieldRP jacobian(grid);
    ddc::for_each(grid, [&](IndexRP const irp) {
        CoordXY coords_xy = mapping(ddc::coordinate(irp));
        coords_x(irp) = ddc::select<RDimX>(coords_xy);
        coords_y(irp) = ddc::select<RDimY>(coords_xy);
        jacobian(irp) = mapping.jacobian(ddc::coordinate(irp));
    });

    double const tau(1e-10);
    double const phi_max(1.);


    VortexMergerEquilibria equilibrium(mapping, grid, builder, spline_evaluator, poisson_solver);
    std::function<double(double const)> const function = [&](double const x) { return x * x; };
    DFieldRP rho_eq(grid);
    equilibrium.set_equilibrium(rho_eq, function, phi_max, tau);


    VortexMergerDensitySolution solution(mapping);
    DFieldRP rho(grid);
    solution.set_initialisation(rho, rho_eq, eps, sigma, x_star_1, y_star_1, x_star_2, y_star_2);


    // Compute phi equilibrium phi_eq from Poisson solver. ***********
    DFieldRP phi_eq(grid);
    Spline2D rho_coef_eq(dom_bsplinesRP);
    builder(rho_coef_eq.span_view(), rho_eq.span_cview());
    PoissonLikeRHSFunction poisson_rhs_eq(rho_coef_eq.span_view(), spline_evaluator);
    poisson_solver(poisson_rhs_eq, coords.span_cview(), phi_eq.span_view());


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
    predcorr_operator(rho.span_view(), dt, iter_nb);


    end_simulation = std::chrono::system_clock::now();
    display_time_difference("Simulation time: ", start_simulation, end_simulation);


    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
    PC_tree_destroy(&conf_gyselalibxx);

    return EXIT_SUCCESS;
}
