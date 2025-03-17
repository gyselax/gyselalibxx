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
#include "spline_interpolator_2d.hpp"
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

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(Nt);

    std::vector<CoordR> r_knots = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> theta_knots
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);

    // Creating mesh & supports:
    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR const interpolation_idx_range_r(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta const interpolation_idx_range_theta(
            SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const grid(interpolation_idx_range_r, interpolation_idx_range_theta);

    // Split the index range of the advection field along RTheta
    const int npoints_theta = IdxRangeTheta(grid).size();
    IdxRangeRTheta const grid_without_Opoint(grid.remove_first(IdxStepRTheta(1, 0)));
    IdxRangeRTheta const Opoint_grid(grid.take_first(IdxStepRTheta(1, npoints_theta)));


    // OPERATORS ======================================================================================
    SplineRThetaBuilder_host const builder_host(grid);
    SplineRThetaBuilder const builder(grid);

    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(r_min);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(r_max);

    SplineRThetaEvaluatorConstBound_host spline_evaluator_extrapol_host(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());
    SplineRThetaEvaluatorConstBound spline_evaluator_extrapol(
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
            spline_evaluator_extrapol_host);
    DiscreteToCartesian const discrete_mapping = discrete_mapping_builder();

    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping);


    // --- Advection operator -------------------------------------------------------------------------
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

    PreallocatableSplineInterpolator2D interpolator(builder, spline_evaluator);

    RK3<FieldMemRTheta<CoordRTheta>,
        DVectorFieldMemRTheta<X, Y>,
        Kokkos::DefaultExecutionSpace> const time_stepper(grid);
    SplinePolarFootFinder find_feet(
            time_stepper,
            to_physical_mapping,
            to_physical_mapping,
            builder,
            spline_evaluator_extrapol);

    BslAdvectionRTheta advection_operator(interpolator, find_feet, to_physical_mapping);

    // --- Advection field finder ---------------------------------------------------------------------
    AdvectionFieldFinder advection_field_computer(to_physical_mapping);


    // --- Choice of the simulation -------------------------------------------------------------------
#if defined(TRANSLATION)
    AdvectionFieldSimulation simulation
            = get_translation_advection_field_simulation(to_physical_mapping, rmin, rmax);
#elif defined(ROTATION)
    AdvectionFieldSimulation simulation
            = get_rotation_advection_field_simulation(to_physical_mapping, rmin, rmax);
#elif defined(DECENTRED_ROTATION)
    AdvectionFieldSimulation simulation
            = get_decentred_rotation_advection_field_simulation(to_physical_mapping);
#endif

    // ================================================================================================
    // SIMULATION DATA                                                                                 |
    // ================================================================================================

    // --- Time parameters ----------------------------------------------------------------------------
    int const iter_nb = final_T * int(1 / dt);

    // ================================================================================================
    // INITIALISATION                                                                                 |
    // ================================================================================================
    host_t<DFieldMemRTheta> allfdistribu_rtheta_alloc(grid);
    host_t<DFieldMemRTheta> allfdistribu_xy_alloc(grid);

    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_exact_alloc(grid);
    host_t<DVectorFieldMemRTheta<R, Theta>> advection_field_rtheta_alloc(grid_without_Opoint);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_alloc(grid);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_from_rtheta_alloc(grid);
    CoordXY advection_field_xy_centre;

    host_t<DFieldMemRTheta> electrostatic_potential_alloc(grid);

    host_t<DFieldRTheta> allfdistribu_rtheta(allfdistribu_rtheta_alloc);
    host_t<DFieldRTheta> allfdistribu_xy(allfdistribu_xy_alloc);
    host_t<DFieldRTheta> electrostatic_potential(electrostatic_potential_alloc);

    host_t<DVectorFieldRTheta<X, Y>> advection_field_exact(advection_field_exact_alloc);
    host_t<DVectorFieldRTheta<R, Theta>> advection_field_rtheta(advection_field_rtheta_alloc);
    host_t<DVectorFieldRTheta<X, Y>> advection_field_xy(advection_field_xy_alloc);
    host_t<DVectorFieldRTheta<X, Y>> advection_field_xy_from_rtheta(
            advection_field_xy_from_rtheta_alloc);



    // Initialise functions ******************************************
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));
        CoordXY const coord_xy(to_physical_mapping(coord_rtheta));

        allfdistribu_rtheta(irtheta) = simulation.function(coord_rtheta);
        allfdistribu_xy(irtheta) = allfdistribu_rtheta(irtheta);
        electrostatic_potential(irtheta) = simulation.electrostatical_potential(coord_xy, 0);

        CoordXY const evaluated_advection_field = simulation.advection_field(coord_xy, 0);
        ddcHelper::get<X>(advection_field_exact)(irtheta) = CoordX(evaluated_advection_field);
        ddcHelper::get<Y>(advection_field_exact)(irtheta) = CoordY(evaluated_advection_field);
    });


    // Constant advection fields *************************************
    advection_field_computer(
            electrostatic_potential,
            advection_field_rtheta,
            advection_field_xy_centre);
    advection_field_computer(electrostatic_potential, advection_field_xy);


    // Compare advection fields ---
    host_t<DVectorFieldMemRTheta<X, Y>> difference_between_fields_exact_and_xy(grid);
    // > Compare the advection field computed on XY to the exact advection field
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        ddcHelper::get<X>(difference_between_fields_exact_and_xy)(irtheta)
                = ddcHelper::get<X>(advection_field_exact)(irtheta)
                  - ddcHelper::get<X>(advection_field_xy)(irtheta);
        ddcHelper::get<Y>(difference_between_fields_exact_and_xy)(irtheta)
                = ddcHelper::get<Y>(advection_field_exact)(irtheta)
                  - ddcHelper::get<Y>(advection_field_xy)(irtheta);
    });


    // > Compare the advection field computed on RTheta to the advection field computed on XY
    host_t<DVectorFieldMemRTheta<X, Y>> difference_between_fields_xy_and_rtheta(grid);

    ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irtheta) {
        CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

        std::array<std::array<double, 2>, 2> J; 
        to_physical_mapping.jacobian_matrix(coord_rtheta, J);

        // computation made in BslAdvectionRTheta operator:
        ddcHelper::get<X>(advection_field_xy_from_rtheta)(irtheta)
                = ddcHelper::get<R>(advection_field_rtheta)(irtheta) * J[0][0]
                  + ddcHelper::get<Theta>(advection_field_rtheta)(irtheta) * J[0][1];
        ddcHelper::get<Y>(advection_field_xy_from_rtheta)(irtheta)
                = ddcHelper::get<R>(advection_field_rtheta)(irtheta) * J[1][0]
                  + ddcHelper::get<Theta>(advection_field_rtheta)(irtheta) * J[1][1];

        // compare
        ddcHelper::get<X>(difference_between_fields_xy_and_rtheta)(irtheta)
                = ddcHelper::get<X>(advection_field_xy_from_rtheta)(irtheta)
                  - ddcHelper::get<X>(advection_field_xy)(irtheta);
        ddcHelper::get<Y>(difference_between_fields_xy_and_rtheta)(irtheta)
                = ddcHelper::get<Y>(advection_field_xy_from_rtheta)(irtheta)
                  - ddcHelper::get<Y>(advection_field_xy)(irtheta);
    });

    ddc::for_each(Opoint_grid, [&](IdxRTheta const irtheta) {
        // computation made in BslAdvectionRTheta operator:
        ddcHelper::get<X>(advection_field_xy_from_rtheta)(irtheta)
                = CoordX(advection_field_xy_centre);
        ddcHelper::get<Y>(advection_field_xy_from_rtheta)(irtheta)
                = CoordY(advection_field_xy_centre);

        // compare
        ddcHelper::get<X>(difference_between_fields_xy_and_rtheta)(irtheta)
                = ddcHelper::get<X>(advection_field_xy_from_rtheta)(irtheta)
                  - ddcHelper::get<X>(advection_field_xy)(irtheta);
        ddcHelper::get<Y>(difference_between_fields_xy_and_rtheta)(irtheta)
                = ddcHelper::get<Y>(advection_field_xy_from_rtheta)(irtheta)
                  - ddcHelper::get<Y>(advection_field_xy)(irtheta);
    });

    // --- Check the difference on advection fields  --------------------------------------------------
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        EXPECT_LE(abs(ddcHelper::get<X>(difference_between_fields_exact_and_xy)(irtheta)), 1e-5);
        EXPECT_LE(abs(ddcHelper::get<Y>(difference_between_fields_exact_and_xy)(irtheta)), 1e-5);

        EXPECT_LE(abs(ddcHelper::get<X>(difference_between_fields_xy_and_rtheta)(irtheta)), 1e-13);
        EXPECT_LE(abs(ddcHelper::get<Y>(difference_between_fields_xy_and_rtheta)(irtheta)), 1e-13);
    });

    auto allfdistribu_xy_device
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), allfdistribu_xy);
    auto advection_field_xy_device = ddcHelper::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), advection_field_xy);

    // ================================================================================================
    // SIMULATION                                                                                     |
    // ================================================================================================
    for (int iter(0); iter < iter_nb; ++iter) {
        advection_operator(
                allfdistribu_rtheta,
                advection_field_rtheta,
                advection_field_xy_centre,
                dt);
        advection_operator(
                get_field(allfdistribu_xy_device),
                get_const_field(advection_field_xy_device),
                dt);

        ddc::parallel_deepcopy(allfdistribu_xy, get_const_field(allfdistribu_xy_device));

        // Check the advected functions ---
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            EXPECT_NEAR(allfdistribu_rtheta(irtheta), allfdistribu_xy(irtheta), 5e-13);
        });
    }

    end_simulation = std::chrono::system_clock::now();
    display_time_difference("Simulation time: ", start_simulation, end_simulation);
}
