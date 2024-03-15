// SPDX-License-Identifier: MIT

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

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
#include "diocotron_initialization_equilibrium.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "poisson_rhs_function.hpp"
#include "polarpoissonsolver.hpp"
#include "quadrature.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "spline_quadrature.hpp"
#include "trapezoid_quadrature.hpp"
#include "vlasovpoissonsolver.hpp"



namespace {
using PoissonSolver = PolarSplineFEMPoissonSolver;
using DiscreteMapping = DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>;
using Mapping = CircularToCartesian<RDimX, RDimY, RDimR, RDimP>;

using Evaluator = SplineEvaluator2D<BSplinesR, BSplinesP>;
using Builder = SplineBuilder2D<SplineRBuilder, SplinePBuilder>;
using Interpolator = SplineInterpolatorRP;

namespace fs = std::filesystem;

} // end namespace



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

    double const W1(PCpp_double(conf_gyselalibxx, ".Mesh.r_min"));
    double const R1(PCpp_double(conf_gyselalibxx, ".Mesh.r_minus"));
    double const R2(PCpp_double(conf_gyselalibxx, ".Mesh.r_plus"));
    double const W2(PCpp_double(conf_gyselalibxx, ".Mesh.r_max"));

    CoordR const r_min(W1);
    CoordR const r_max(W2);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordP> p_knots(p_size + 1);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    r_knots[0] = r_min;
    for (int i(1); i < r_size; ++i) {
        r_knots[i] = r_min + i * dr;
    }
    r_knots[r_size] = r_max;

    p_knots[p_size] = p_min;
    for (int i(1); i < p_size; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }
    p_knots[p_size] = p_max;


    // Creating mesh & supports:
    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<IDimR>(SplineInterpPointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(SplineInterpPointsP::get_sampling());

    IDomainR const interpolation_domain_R(SplineInterpPointsR::get_domain());
    IDomainP const interpolation_domain_P(SplineInterpPointsP::get_domain());
    IDomainRP const grid(interpolation_domain_R, interpolation_domain_P);

    FieldRP<CoordRP> coords(grid);
    ddc::for_each(grid, [&](IndexRP const irp) { coords(irp) = ddc::coordinate(irp); });


    // OPERATORS ======================================================================================
    SplineRBuilder const r_builder(interpolation_domain_R);
    SplinePBuilder const p_builder(interpolation_domain_P);
    SplineRPBuilder const builder(grid);

    // --- Define the mapping. ------------------------------------------------------------------------
    ConstantExtrapolationBoundaryValue2D<BSplinesR, BSplinesP, RDimR> boundary_condition_r_left(
            r_min);
    ConstantExtrapolationBoundaryValue2D<BSplinesR, BSplinesP, RDimR> boundary_condition_r_right(
            r_max);

    SplineRPEvaluator spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>);

    const Mapping mapping;
    DiscreteMapping const discrete_mapping
            = DiscreteMapping::analytical_to_discrete(mapping, builder, spline_evaluator_extrapol);

    ddc::init_discrete_space<PolarBSplinesRP>(discrete_mapping, r_builder, p_builder);

    BSDomainRP const dom_bsplinesRP = builder.spline_domain();


    // --- Time integration method --------------------------------------------------------------------
#if defined(EULER_METHOD)
    Euler<FieldRP<CoordRP>, VectorDFieldRP<RDimX, RDimY>> const time_stepper(grid);

#elif defined(CRANK_NICOLSON_METHOD)
    double const epsilon_CN = 1e-8;
    CrankNicolson<FieldRP<CoordRP>, VectorDFieldRP<RDimX, RDimY>> const
            time_stepper(grid, 20, epsilon_CN);

#elif defined(RK3_METHOD)
    RK3<FieldRP<CoordRP>, VectorDFieldRP<RDimX, RDimY>> const time_stepper(grid);

#elif defined(RK4_METHOD)
    RK4<FieldRP<CoordRP>, VectorDFieldRP<RDimX, RDimY>> const time_stepper(grid);

#endif


    // --- Advection operator -------------------------------------------------------------------------
    Evaluator spline_evaluator(
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>);

    PreallocatableSplineInterpolatorRP interpolator(builder, spline_evaluator);

    AdvectionPhysicalDomain advection_domain(mapping);

    SplineFootFinder
            find_feet(time_stepper, advection_domain, grid, builder, spline_evaluator_extrapol);

    BslAdvectionRP advection_operator(interpolator, find_feet);



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

    builder(coeff_alpha_spline, coeff_alpha);
    builder(coeff_beta_spline, coeff_beta);

    PoissonSolver poisson_solver(coeff_alpha_spline, coeff_beta_spline, discrete_mapping);

    // --- Predictor corrector operator ---------------------------------------------------------------
#if defined(PREDCORR)
    BslPredCorrRP predcorr_operator(
            advection_domain,
            mapping,
            advection_operator,
            builder,
            spline_evaluator,
            poisson_solver);
#elif defined(EXPLICIT_PREDCORR)
    BslExplicitPredCorrRP predcorr_operator(
            advection_domain,
            mapping,
            advection_operator,
            grid,
            builder,
            spline_evaluator,
            poisson_solver,
            spline_evaluator_extrapol);
#elif defined(IMPLICIT_PREDCORR)
    BslImplicitPredCorrRP predcorr_operator(
            advection_domain,
            mapping,
            advection_operator,
            grid,
            builder,
            spline_evaluator,
            poisson_solver,
            spline_evaluator_extrapol);
#endif

    // ================================================================================================
    // SIMULATION DATA                                                                                 |
    // ================================================================================================
    double const Q(PCpp_double(
            conf_gyselalibxx,
            ".Perturbation.charge_Q")); // no charge carried by the inner conductor r = W1.
    int const l(PCpp_int(conf_gyselalibxx, ".Perturbation.l_mode"));
    double const eps(PCpp_double(conf_gyselalibxx, ".Perturbation.eps"));
    DiocotronDensitySolution exact_rho(W1, R1, R2, W2, Q, l, eps);

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

    ddc::expose_to_pdi("slope", exact_rho.get_slope());



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



    DFieldRP rho(grid);
    DFieldRP rho_eq(grid);

    // Initialize rho and rho equilibrium ****************************
    for_each(grid, [&](IndexRP const irp) {
        rho(irp) = exact_rho.initialisation(coords(irp));
        rho_eq(irp) = exact_rho.equilibrium(coords(irp));
    });

    // Compute phi equilibrium phi_eq from Poisson solver. ***********
    DFieldRP phi_eq(grid);
    Spline2D rho_coef_eq(dom_bsplinesRP);
    builder(rho_coef_eq, rho_eq);
    PoissonRHSFunction poisson_rhs_eq(rho_coef_eq, spline_evaluator);
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
