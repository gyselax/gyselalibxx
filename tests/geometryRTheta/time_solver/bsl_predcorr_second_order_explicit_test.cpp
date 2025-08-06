// SPDX-License-Identifier: MIT
#include <chrono>
// #include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include <gtest/gtest.h>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_polar.hpp"
#include "bsl_predcorr.hpp"
#include "bsl_predcorr_second_order_explicit.hpp"
#include "bsl_predcorr_second_order_implicit.hpp"
#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "diocotron_initialisation_equilibrium.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "l_norm_tools.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
// #include "params.yaml.hpp"
// #include "pdi_out.yml.hpp"
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
// #include "quadrature.hpp"
// #include "spline_quadrature.hpp"



namespace {
using PoissonSolver = PolarSplineFEMPoissonLikeSolver<
        GridR,
        GridTheta,
        PolarBSplinesRTheta,
        SplineRThetaEvaluatorNullBound>;
using DiscreteMappingBuilder
        = DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>;
using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;


using SplineRBuilder = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::HostSpace,
        BSplinesR,
        GridR,
        ddc::BoundCond::GREVILLE, // boundary at r=0
        ddc::BoundCond::GREVILLE, // boundary at r_max
        ddc::SplineSolver::LAPACK>;
using SplineThetaBuilder = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::HostSpace,
        BSplinesTheta,
        GridTheta,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK>;


void RMSE(double& slope, double& offset, std::vector<double> const& x, std::vector<double> const& y)
{
    assert(x.size() == y.size());
    std::size_t size = x.size();

    double sum_x = 0;
    double sum_x2 = 0;
    double sum_y = 0;
    double sum_xy = 0;

    for (int i(0); i <= size; i++) {
        sum_x = sum_x + x[i];
        sum_x2 = sum_x2 + x[i] * x[i];
        sum_y = sum_y + y[i];
        sum_xy = sum_xy + x[i] * y[i];
    }
    /* Calculating a and b */
    offset = (size * sum_xy - sum_x * sum_y) / (size * sum_x2 - sum_x * sum_x);
    slope = (sum_y - offset * sum_x) / size;
}
// namespace fs = std::filesystem;


constexpr char const* const PDI_CFG = R"PDI_CFG()PDI_CFG"; 
} // end namespace



TEST(BslPredCorrSecondOrderExplicitTest, CheckGrowthRate)
{
    // SETUP ==========================================================================================
    // fs::create_directory("output");

    // Get the parameters of the mesh_rtheta from the grid_size.yaml. ----------------------------------------
    // PC_tree_t conf_gyselalibxx = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    // Kokkos::ScopeGuard kokkos_scope(argc, argv);
    // ddc::ScopeGuard ddc_scope(argc, argv);

    // std::chrono::time_point<std::chrono::system_clock> start_simulation;
    // std::chrono::time_point<std::chrono::system_clock> end_simulation;

    // start_simulation = std::chrono::system_clock::now();

    // Build the mesh_rtheta for the space. ------------------------------------------------------------------
    CoordR r_min(0.);
    CoordR r_max(1.);
    IdxStepR n_rcells(128 / 4);

    CoordTheta theta_min(0.);
    CoordTheta theta_max(2*M_PI);
    IdxStepTheta n_thetacells(256 / 4);

    std::vector<CoordR> break_points_r = build_uniform_break_points(r_min, r_max, n_rcells);
    std::vector<CoordTheta> break_points_theta
            = build_uniform_break_points(theta_min, theta_max, n_thetacells);

    ddc::init_discrete_space<BSplinesR>(break_points_r);
    ddc::init_discrete_space<BSplinesTheta>(break_points_theta);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::template get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(
            SplineInterpPointsTheta::template get_sampling<GridTheta>());

    IdxRangeR const mesh_r = (SplineInterpPointsR::template get_domain<GridR>());
    IdxRangeTheta const mesh_theta = (SplineInterpPointsTheta::template get_domain<GridTheta>());
    double const dt(0.1);
    double const final_T(40.0);

    IdxRangeRTheta const mesh_rtheta(mesh_r, mesh_theta);

    host_t<FieldMemRTheta<CoordRTheta>> coords(mesh_rtheta);
    ddc::for_each(mesh_rtheta, [&](IdxRTheta const irtheta) {
        coords(irtheta) = ddc::coordinate(irtheta);
    });


    // OPERATORS ======================================================================================
    SplineRThetaBuilder const builder(mesh_rtheta);
    SplineRThetaBuilder_host const builder_host(mesh_rtheta);

    SplineRBuilder const r_builder(mesh_r);
    // SplineRThetaBuilder_host const builder_host(mesh_rtheta);

    SplineThetaBuilder const theta_builder(mesh_theta);
    // SplineRThetaBuilder_host const builder_host(mesh_rtheta);


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
    DiscreteToCartesian const discrete_mapping = discrete_mapping_builder();


    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping);

    IdxRangeBSRTheta const idx_range_bsplinesRTheta = get_spline_idx_range(builder);


    // --- Time integration method --------------------------------------------------------------------
    EulerBuilder const time_stepper;

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

    PreallocatableSplineInterpolator2D interpolator(builder, spline_evaluator, mesh_rtheta);

    SplinePolarFootFinder find_feet(
            mesh_rtheta,
            time_stepper,
            to_physical_mapping,
            to_physical_mapping,
            builder,
            spline_evaluator_extrapol);

    BslAdvectionPolar advection_operator(interpolator, find_feet, to_physical_mapping);



    // --- Poisson solver -----------------------------------------------------------------------------
    // Coefficients alpha and beta of the Poisson equation:
    DFieldMemRTheta coeff_alpha(mesh_rtheta);
    DFieldMemRTheta coeff_beta(mesh_rtheta);
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
    BslExplicitPredCorrRTheta predcorr_operator(
            to_physical_mapping,
            to_physical_mapping,
            advection_operator,
            mesh_rtheta,
            builder_host,
            poisson_solver,
            spline_evaluator_extrapol_host);

    // ================================================================================================
    // SIMULATION DATA                                                                                 |
    // ================================================================================================
    double const Q(0); // no charge carried by the inner conductor r = W1.
    int const l(9);
    double const eps(0.0001);
    CoordR const R1(0.45);
    CoordR const R2(0.50);
    DiocotronDensitySolution exact_rho(
            ddc::coordinate(mesh_r.front()),
            R1,
            R2,
            ddc::coordinate(mesh_r.back()),
            Q,
            l,
            eps);

    const double exact_growth_rate = exact_rho.get_slope();

    // --- Time parameters ----------------------------------------------------------------------------
    int const iter_nb = final_T * int(1 / dt);

    // --- save simulation data
    // ddc::expose_to_pdi("r_size", ddc::discrete_space<BSplinesR>().ncells());
    // ddc::expose_to_pdi("theta_size", ddc::discrete_space<BSplinesTheta>().ncells());

    // expose_mesh_to_pdi("r_coords", mesh_r);
    // expose_mesh_to_pdi("theta_coords", mesh_theta);

    // ddc::expose_to_pdi("delta_t", dt);
    // ddc::expose_to_pdi("final_T", final_T);
    // ddc::expose_to_pdi("time_step_diag", PCpp_int(conf_gyselalibxx, ".Output.time_step_diag"));

    // ddc::expose_to_pdi("slope", exact_rho.get_slope());

    // --- Quadrature to compute the L2-norm
    host_t<DFieldMemRTheta> const quadrature_coeffs = compute_coeffs_on_mapping(
            Kokkos::DefaultHostExecutionSpace(),
            to_physical_mapping,
            spline_quadrature_coefficients<
                    Kokkos::DefaultHostExecutionSpace>(mesh_rtheta, r_builder, theta_builder));
    host_t<Quadrature<IdxRangeRTheta>> quadrature(get_const_field(quadrature_coeffs));

    // ================================================================================================
    // INITIALISATION                                                                                 |
    // ================================================================================================
    // Cartesian coordinates and jacobian ****************************
    host_t<FieldMemRTheta<CoordX>> coords_x(mesh_rtheta);
    host_t<FieldMemRTheta<CoordY>> coords_y(mesh_rtheta);
    host_t<DFieldMemRTheta> jacobian(mesh_rtheta);
    ddc::for_each(mesh_rtheta, [&](IdxRTheta const irtheta) {
        CoordXY coords_xy = to_physical_mapping(ddc::coordinate(irtheta));
        coords_x(irtheta) = ddc::select<X>(coords_xy);
        coords_y(irtheta) = ddc::select<Y>(coords_xy);
        jacobian(irtheta) = to_physical_mapping.jacobian(ddc::coordinate(irtheta));
    });

    host_t<DFieldMemRTheta> rho_alloc(mesh_rtheta);
    host_t<DFieldMemRTheta> rho_eq_alloc(mesh_rtheta);

    auto rho_alloc_host = ddc::create_mirror_view(Kokkos::DefaultHostExecutionSpace(), get_field(rho_alloc));

    // Initialise rho and rho equilibrium ****************************
    ddc::for_each(mesh_rtheta, [&](IdxRTheta const irtheta) {
        rho_alloc(irtheta) = exact_rho.initialisation(coords(irtheta));
        rho_eq_alloc(irtheta) = exact_rho.equilibrium(coords(irtheta));
    });

    // Compute phi equilibrium phi_eq from Poisson solver. ***********
    DFieldMemRTheta phi_eq(mesh_rtheta);
    host_t<DFieldMemRTheta> phi_eq_host(mesh_rtheta);
    host_t<Spline2DMem> rho_coef_eq(idx_range_bsplinesRTheta);
    builder_host(get_field(rho_coef_eq), get_const_field(rho_eq_alloc));
    PoissonLikeRHSFunction poisson_rhs_eq(get_const_field(rho_coef_eq), spline_evaluator_host);
    poisson_solver(poisson_rhs_eq, get_field(phi_eq));
    ddc::parallel_deepcopy(phi_eq, phi_eq_host);

    // // --- Save initial data --------------------------------------------------------------------------
    // ddc::PdiEvent("initialisation")
    //         .with("x_coords", coords_x)
    //         .with("y_coords", coords_y)
    //         .with("jacobian", jacobian)
    //         .with("density_eq", rho_eq_alloc)
    //         .with("electrical_potential_eq", phi_eq_host);


    // ================================================================================================
    // SIMULATION                                                                                     |
    // ================================================================================================
    std::vector<double> L2_norms;
    std::vector<double> time_steps;
    for (int iter(0); iter < iter_nb; iter++) {
        // compute the density at the next time step
        predcorr_operator(get_field(rho_alloc), dt, iter);
        // compute the L2 norm
        start_time = std::chrono::system_clock::now();
        ddc::parallel_deepcopy(get_field(rho_alloc_host), get_field(rho_alloc));
        L2_norms.push_back(
                norm_L2(Kokkos::DefaultHostExecutionSpace(), quadrature, get_field(rho_alloc_host)));
        time_steps.push_back(dt * iter);
        end_time = std::chrono::system_clock::now();
        display_time_difference(">>> Computation L2_norm time: ", start_time, end_time);
    }


    // end_simulation = std::chrono::system_clock::now();
    // display_time_difference("Simulation time: ", start_simulation, end_simulation);

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
    // PC_tree_destroy(&conf_gyselalibxx);

    // Check the growth rate.
    double offset;
    double slope;
    RMSE(slope, offset, time_steps, L2_norms);


    std::cout << "slope = " << slope 
    << "     exact_growth_rate = " << exact_growth_rate 
    << "     diff = " << exact_growth_rate -slope 
    << std::endl; 
    
    EXPECT_NEAR(slope, exact_growth_rate, 1e-13); 


    // return EXIT_SUCCESS;
}
