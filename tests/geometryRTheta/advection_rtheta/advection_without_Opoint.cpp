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

#include "../advection_field_rtheta/test_cases_adv_field.hpp"

#include "bsl_advection_polar.hpp"
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
#include "inverse_jacobian_matrix.hpp"
#include "l_norm_tools.hpp"
#include "mesh_builder.hpp"
#include "quadrature.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_interpolator_2d.hpp"
#include "spline_polar_foot_finder.hpp"
#include "spline_quadrature.hpp"
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



TEST(AdvectionWithoutOpointComputation, TestAdvectionFieldFinder)
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

    double const rmin(0.1);
    double const rmax(1);

    CoordR const r_min(rmin);
    CoordR const r_max(rmax);
    IdxStepR const r_ncells(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(Nt);

    std::vector<CoordR> r_break_points = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> theta_break_points
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);

    // Creating mesh & supports:
    ddc::init_discrete_space<BSplinesR>(r_break_points);
    ddc::init_discrete_space<BSplinesTheta>(theta_break_points);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR const interpolation_idx_range_r(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta const interpolation_idx_range_theta(
            SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const grid(interpolation_idx_range_r, interpolation_idx_range_theta);

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

    PreallocatableSplineInterpolator2D interpolator(builder, spline_evaluator, grid);

    RK3Builder const time_stepper;
    SplinePolarFootFinder find_feet(
            grid,
            time_stepper,
            to_physical_mapping,
            to_physical_mapping,
            builder,
            spline_evaluator_extrapol);

    BslAdvectionPolar advection_operator(interpolator, find_feet, to_physical_mapping);

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
    host_t<DFieldMemRTheta> density_rtheta_averaged_alloc(grid);
    host_t<DFieldMemRTheta> density_xy_alloc(grid);

    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_exact_alloc(grid);
    host_t<DVectorFieldMemRTheta<R, Theta>> advection_field_rtheta_alloc(grid);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_alloc(grid);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_from_rtheta_alloc(grid);

    host_t<DFieldMemRTheta> electrostatic_potential_alloc(grid);

    host_t<DFieldRTheta> density_rtheta_averaged(density_rtheta_averaged_alloc);
    host_t<DFieldRTheta> density_xy(density_xy_alloc);
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

        density_rtheta_averaged(irtheta) = simulation.function(coord_rtheta);
        density_xy(irtheta) = density_rtheta_averaged(irtheta);
        electrostatic_potential(irtheta) = simulation.electrostatical_potential(coord_xy, 0);

        ddcHelper::assign_vector_field_element(
                advection_field_exact,
                irtheta,
                simulation.advection_field(coord_xy, 0));
    });


    // Constant advection fields *************************************
    advection_field_computer(electrostatic_potential, advection_field_xy);

    InverseJacobianMatrix inv_jacobian_matrix(to_physical_mapping);
    ddc::for_each(grid, [&](Idx<GridR, GridTheta> const idx) {
        Coord<R, Theta> const coord_rtheta(ddc::coordinate(idx));
        Tensor inv_J = inv_jacobian_matrix(coord_rtheta);

        ddcHelper::assign_vector_field_element(
                advection_field_rtheta,
                idx,
                tensor_mul(index<'i', 'j'>(inv_J), index<'j'>(advection_field_xy(idx))));
    });

    auto density_xy_device
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), density_xy);
    auto advection_field_xy_device = ddcHelper::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), advection_field_xy);


    auto density_rtheta_averaged_device = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), density_rtheta_averaged);
    auto advection_field_rtheta_device = ddcHelper::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), advection_field_rtheta);

    // ================================================================================================
    // SIMULATION                                                                                     |
    // ================================================================================================
    for (int iter(0); iter < iter_nb; ++iter) {
        // --- operator() 2: compute a value for the O-point from the other values.
        advection_operator(
                get_field(density_rtheta_averaged_device),
                get_const_field(advection_field_rtheta_device),
                dt);
        // --- operator() 3: directly give the advection field on (x,y). No extra computations.
        advection_operator(
                get_field(density_xy_device),
                get_const_field(advection_field_xy_device),
                dt);

        ddc::parallel_deepcopy(density_xy, get_const_field(density_xy_device));
        ddc::parallel_deepcopy(
                density_rtheta_averaged,
                get_const_field(density_rtheta_averaged_device));

        // Check the advected functions ---
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            EXPECT_NEAR(density_rtheta_averaged(irtheta), density_xy(irtheta), 5e-7);
        });
    }

    end_simulation = std::chrono::system_clock::now();
    display_time_difference("Simulation time: ", start_simulation, end_simulation);
}
