/**
 * Test the computation of the advection field passing by polar axis. 
 * Also test the advection with an advection field along the polar axis given as input. 
*/

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_rtheta.hpp"
#include "bsl_predcorr.hpp"
#include "bsl_predcorr_second_order_explicit.hpp"
#include "bsl_predcorr_second_order_implicit.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "l_norm_tools.hpp"
#include "mesh_builder.hpp"
#include "quadrature.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "spline_polar_foot_finder.hpp"
#include "spline_quadrature.hpp"
#include "test_cases_adv_field.hpp"
#include "trapezoid_quadrature.hpp"


namespace {
using DiscreteMappingBuilder = DiscreteToCartesianBuilder<
        X,
        Y,
        SplineRThetaBuilder_host,
        SplineRThetaEvaluatorConstBound_host>;
using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;

namespace fs = std::filesystem;

} // end namespace



TEST(AdvectionFieldRThetaComputation, TestAdvectionFieldFinder)
{
    // SETUP ==========================================================================================
    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;

    start_simulation = std::chrono::system_clock::now();
    // Build the grid for the space. ------------------------------------------------------------------
    int const Nr(20);
    int const Nt(40);
    double const dt(0.1);
    double const final_T(0.8);

    double const rmin(0);
    double const rmax(1);

    CoordR const r_min(rmin);
    CoordR const r_max(rmax);
    IdxStepR const r_ncells(Nr);

    CoordTheta const p_min(0.0);
    CoordTheta const p_max(2.0 * M_PI);
    IdxStepTheta const p_ncells(Nt);

    std::vector<CoordR> r_knots = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> p_knots = build_uniform_break_points(p_min, p_max, p_ncells);

    // Creating mesh & supports:
    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(p_knots);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR const interpolation_idx_range_R(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta const interpolation_idx_range_P(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const grid(interpolation_idx_range_R, interpolation_idx_range_P);

    // Split the index range of the advection field along RTheta
    const int npoints_p = IdxRangeTheta(grid).size();
    IdxRangeRTheta const grid_without_Opoint(grid.remove_first(IdxStepRTheta(1, 0)));
    IdxRangeRTheta const Opoint_grid(grid.take_first(IdxStepRTheta(1, npoints_p)));


    // OPERATORS ======================================================================================
    SplineRThetaBuilder_host const builder_host(grid);
    SplineRThetaBuilder const builder(grid);

    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(r_min);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(r_max);

    SplineRThetaEvaluatorConstBound_host spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());


    ddc::NullExtrapolationRule r_extrapolation_rule;
    PolarSplineEvaluator<PolarBSplinesRTheta, ddc::NullExtrapolationRule> polar_spline_evaluator(
            r_extrapolation_rule);

    // --- Define the to_physical_mapping. ------------------------------------------------------------------------
    const LogicalToPhysicalMapping to_physical_mapping;
    DiscreteMappingBuilder const discrete_mapping_builder(
            Kokkos::DefaultHostExecutionSpace(),
            to_physical_mapping,
            builder_host,
            spline_evaluator_extrapol);
    DiscreteToCartesian const discrete_mapping = discrete_mapping_builder();

    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping);


    // --- Advection operator -------------------------------------------------------------------------
    ddc::PeriodicExtrapolationRule<Theta> p_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    PreallocatableSplineInterpolatorRTheta interpolator(builder, spline_evaluator);

    RK3<host_t<FieldMemRTheta<CoordRTheta>>,
        host_t<DVectorFieldMemRTheta<X, Y>>,
        Kokkos::DefaultHostExecutionSpace> const time_stepper(grid);
    SplinePolarFootFinder find_feet(
            time_stepper,
            to_physical_mapping,
            to_physical_mapping,
            builder_host,
            spline_evaluator_extrapol);

    BslAdvectionRTheta advection_operator(interpolator, find_feet, to_physical_mapping);

    // --- Advection field finder ---------------------------------------------------------------------
    AdvectionFieldFinder advection_field_computer(to_physical_mapping);


    // --- Choice of the simulation -------------------------------------------------------------------
#if defined(TRANSLATION)
    TranslationAdvectionFieldSimulation simulation(to_physical_mapping, rmin, rmax);
#elif defined(ROTATION)
    RotationAdvectionFieldSimulation simulation(to_physical_mapping, rmin, rmax);
#elif defined(DECENTRED_ROTATION)
    DecentredRotationAdvectionFieldSimulation simulation(to_physical_mapping);
#endif

    // ================================================================================================
    // SIMULATION DATA                                                                                 |
    // ================================================================================================

    // --- Time parameters ----------------------------------------------------------------------------
    int const iter_nb = final_T * int(1 / dt);

    // ================================================================================================
    // INITIALISATION                                                                                 |
    // ================================================================================================
    host_t<DFieldMemRTheta> allfdistribu_rp_alloc(grid);
    host_t<DFieldMemRTheta> allfdistribu_xy_alloc(grid);

    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_exact_alloc(grid);
    host_t<DVectorFieldMemRTheta<R, Theta>> advection_field_rtheta_alloc(grid_without_Opoint);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_alloc(grid);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_from_rp_alloc(grid);
    CoordXY advection_field_xy_center;

    host_t<DFieldMemRTheta> electrostatic_potential_alloc(grid);

    host_t<DFieldRTheta> allfdistribu_rp(allfdistribu_rp_alloc);
    host_t<DFieldRTheta> allfdistribu_xy(allfdistribu_xy_alloc);
    host_t<DFieldRTheta> electrostatic_potential(electrostatic_potential_alloc);

    host_t<DVectorFieldRTheta<X, Y>> advection_field_exact(advection_field_exact_alloc);
    host_t<DVectorFieldRTheta<R, Theta>> advection_field_rtheta(advection_field_rtheta_alloc);
    host_t<DVectorFieldRTheta<X, Y>> advection_field_xy(advection_field_xy_alloc);
    host_t<DVectorFieldRTheta<X, Y>> advection_field_xy_from_rp(advection_field_xy_from_rp_alloc);



    // Initialize functions ******************************************
    auto function = simulation.get_test_function();
    auto phi_function = simulation.get_electrostatique_potential();
    auto advection_field = simulation.get_advection_field();
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        CoordRTheta const coord_rp(ddc::coordinate(irp));
        CoordXY const coord_xy(to_physical_mapping(coord_rp));

        allfdistribu_rp(irp) = function(coord_rp);
        allfdistribu_xy(irp) = allfdistribu_rp(irp);
        electrostatic_potential(irp) = phi_function(coord_xy, 0);

        CoordXY const evaluated_advection_field = advection_field(coord_xy, 0);
        ddcHelper::get<X>(advection_field_exact)(irp) = CoordX(evaluated_advection_field);
        ddcHelper::get<Y>(advection_field_exact)(irp) = CoordY(evaluated_advection_field);
    });


    // Constant advection fields *************************************
    advection_field_computer(
            electrostatic_potential,
            advection_field_rtheta,
            advection_field_xy_center);
    advection_field_computer(electrostatic_potential, advection_field_xy);


    // Compare advection fields ---
    host_t<DVectorFieldMemRTheta<X, Y>> difference_between_fields_exact_and_xy(grid);
    // > Compare the advection field computed on XY to the exact advection field
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        ddcHelper::get<X>(difference_between_fields_exact_and_xy)(irp)
                = ddcHelper::get<X>(advection_field_exact)(irp)
                  - ddcHelper::get<X>(advection_field_xy)(irp);
        ddcHelper::get<Y>(difference_between_fields_exact_and_xy)(irp)
                = ddcHelper::get<Y>(advection_field_exact)(irp)
                  - ddcHelper::get<Y>(advection_field_xy)(irp);
    });


    // > Compare the advection field computed on RTheta to the advection field computed on XY
    host_t<DVectorFieldMemRTheta<X, Y>> difference_between_fields_xy_and_rp(grid);

    MetricTensor<LogicalToPhysicalMapping, CoordRTheta> metric_tensor(to_physical_mapping);
    ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irp) {
        CoordRTheta const coord_rp(ddc::coordinate(irp));

        std::array<std::array<double, 2>, 2> inv_J; // inverse Jacobian matrix
        to_physical_mapping.inv_jacobian_matrix(coord_rp, inv_J);
        double const jacobian = to_physical_mapping.jacobian(coord_rp);

        // computation made in BslAdvectionRTheta operator:
        ddcHelper::get<X>(advection_field_xy_from_rp)(irp)
                = ddcHelper::get<R>(advection_field_rtheta)(irp) * inv_J[0][0] * jacobian 
                  + ddcHelper::get<Theta>(advection_field_rtheta)(irp) * inv_J[1][0] * jacobian; 
        ddcHelper::get<Y>(advection_field_xy_from_rp)(irp)
                = ddcHelper::get<R>(advection_field_rtheta)(irp) * inv_J[0][1] *jacobian
                  + ddcHelper::get<Theta>(advection_field_rtheta)(irp) * inv_J[1][1] *jacobian; 

        // compare
        ddcHelper::get<X>(difference_between_fields_xy_and_rp)(irp)
                = ddcHelper::get<X>(advection_field_xy_from_rp)(irp)
                  - ddcHelper::get<X>(advection_field_xy)(irp);
        ddcHelper::get<Y>(difference_between_fields_xy_and_rp)(irp)
                = ddcHelper::get<Y>(advection_field_xy_from_rp)(irp)
                  - ddcHelper::get<Y>(advection_field_xy)(irp);
    });

    ddc::for_each(Opoint_grid, [&](IdxRTheta const irp) {
        // computation made in BslAdvectionRTheta operator:
        ddcHelper::get<X>(advection_field_xy_from_rp)(irp) = CoordX(advection_field_xy_center);
        ddcHelper::get<Y>(advection_field_xy_from_rp)(irp) = CoordY(advection_field_xy_center);

        // compare
        ddcHelper::get<X>(difference_between_fields_xy_and_rp)(irp)
                = ddcHelper::get<X>(advection_field_xy_from_rp)(irp)
                  - ddcHelper::get<X>(advection_field_xy)(irp);
        ddcHelper::get<Y>(difference_between_fields_xy_and_rp)(irp)
                = ddcHelper::get<Y>(advection_field_xy_from_rp)(irp)
                  - ddcHelper::get<Y>(advection_field_xy)(irp);
    });

    // --- Check the difference on advection fields  --------------------------------------------------
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        EXPECT_LE(abs(ddcHelper::get<X>(difference_between_fields_exact_and_xy)(irp)), 1e-5);
        EXPECT_LE(abs(ddcHelper::get<Y>(difference_between_fields_exact_and_xy)(irp)), 1e-5);

        EXPECT_LE(abs(ddcHelper::get<X>(difference_between_fields_xy_and_rp)(irp)), 1e-13);
        EXPECT_LE(abs(ddcHelper::get<Y>(difference_between_fields_xy_and_rp)(irp)), 1e-13);
    });


    // ================================================================================================
    // SIMULATION                                                                                     |
    // ================================================================================================
    for (int iter(0); iter < iter_nb; ++iter) {
        advection_operator(allfdistribu_rp, advection_field_rtheta, advection_field_xy_center, dt);
        advection_operator(allfdistribu_xy, advection_field_xy, dt);

        // Check the advected functions ---
        ddc::for_each(grid, [&](IdxRTheta const irp) {
            EXPECT_NEAR(allfdistribu_rp(irp), allfdistribu_xy(irp), 1e-13);
        });
    }

    end_simulation = std::chrono::system_clock::now();
    display_time_difference("Simulation time: ", start_simulation, end_simulation);
}
