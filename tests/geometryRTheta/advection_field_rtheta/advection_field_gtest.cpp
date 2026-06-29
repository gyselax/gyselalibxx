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

#include "bsl_advection_polar.hpp"
#include "bsl_predcorr.hpp"
#include "bsl_predcorr_second_order_explicit.hpp"
#include "bsl_predcorr_second_order_implicit.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_poloidal_cs_spline_mapping.hpp"
#include "discrete_poloidal_cs_spline_mapping_builder.hpp"
#include "euler.hpp"
#include "geometry_r_theta.hpp"
#include "l_norm_tools.hpp"
#include "mesh_builder.hpp"
#include "polar_foot_finder.hpp"
#include "quadrature.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "simulation_utils_tools.hpp"
#include "spline_definitions_r_theta.hpp"
#include "spline_quadrature.hpp"
#include "test_cases_adv_field.hpp"
#include "trapezoid_quadrature.hpp"


namespace {
using DiscreteMappingBuilder = DiscretePoloidalCSSplineMappingBuilder<
        X,
        Y,
        SplineRThetaBuilder_host,
        SplineRThetaEvaluatorConstBound_host>;
using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
using PhysicalToLogicalMapping = CartesianToCircular<X, Y, R, Theta>;

enum class SimulationType { TRANSLATION, ROTATION, DECENTRED_ROTATION };

template <SimulationType sim>
auto get_simulation(
        LogicalToPhysicalMapping to_physical_mapping,
        CoordR const rmin,
        CoordR const rmax)
{
    if constexpr (sim == SimulationType::TRANSLATION) {
        return get_translation_advection_field_simulation(to_physical_mapping, rmin, rmax);
    } else if constexpr (sim == SimulationType::ROTATION) {
        return get_rotation_advection_field_simulation(to_physical_mapping, rmin, rmax);
    } else if constexpr (sim == SimulationType::DECENTRED_ROTATION) {
        return get_decentred_rotation_advection_field_simulation(to_physical_mapping);
    }
}

namespace fs = std::filesystem;

template <class T>
struct AdvectionFieldRThetaComputationFixture;

template <
        SimulationType simulation_choice_,
        FootFindingSpace FFSpace_,
        AdvectionFieldSpace AFSpace_,
        int dt_inv>
struct AdvectionFieldRThetaComputationFixture<std::tuple<
        std::integral_constant<SimulationType, simulation_choice_>,
        std::integral_constant<FootFindingSpace, FFSpace_>,
        std::integral_constant<AdvectionFieldSpace, AFSpace_>,
        std::integral_constant<int, dt_inv>>> : public testing::Test
{
    static constexpr SimulationType simulation_choice = simulation_choice_;
    static constexpr FootFindingSpace FFSpace = FFSpace_;
    static constexpr AdvectionFieldSpace AFSpace = AFSpace_;
    static constexpr double dt = 1.0 / dt_inv;
};

using Cases = ::testing::Types<
        std::tuple<
                std::integral_constant<SimulationType, SimulationType::TRANSLATION>,
                std::integral_constant<FootFindingSpace, FootFindingSpace::PSEUDO_PHYSICAL>,
                std::integral_constant<AdvectionFieldSpace, AdvectionFieldSpace::PHYSICAL>,
                std::integral_constant<int, 100>>,
        std::tuple<
                std::integral_constant<SimulationType, SimulationType::ROTATION>,
                std::integral_constant<FootFindingSpace, FootFindingSpace::PSEUDO_PHYSICAL>,
                std::integral_constant<AdvectionFieldSpace, AdvectionFieldSpace::PHYSICAL>,
                std::integral_constant<int, 100>>,
        std::tuple<
                std::integral_constant<SimulationType, SimulationType::DECENTRED_ROTATION>,
                std::integral_constant<FootFindingSpace, FootFindingSpace::PSEUDO_PHYSICAL>,
                std::integral_constant<AdvectionFieldSpace, AdvectionFieldSpace::PHYSICAL>,
                std::integral_constant<int, 10000>>,
        std::tuple<
                std::integral_constant<SimulationType, SimulationType::ROTATION>,
                std::integral_constant<FootFindingSpace, FootFindingSpace::LOGICAL>,
                std::integral_constant<AdvectionFieldSpace, AdvectionFieldSpace::LOGICAL>,
                std::integral_constant<int, 10>>>;

} // end namespace

TYPED_TEST_SUITE(AdvectionFieldRThetaComputationFixture, Cases);

TYPED_TEST(AdvectionFieldRThetaComputationFixture, TestAdvectionFieldFinder)
{
    // SETUP ==========================================================================================
    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;

    start_simulation = std::chrono::system_clock::now();
    // Build the grid for the space. ------------------------------------------------------------------
    int const Nr(20);
    int const Nt(40);
    double const dt(TestFixture::dt);
    int const iter_nb = 8;

    double const rmin(0);
    double const rmax(1);

    CoordR const r_min(rmin);
    CoordR const r_max(rmax);
    IdxStepR const r_ncells(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(Nt);

    std::vector<CoordR> r_break_points = build_uniform_break_points(CoordR(0), r_max, r_ncells);
    std::vector<CoordTheta> theta_break_points
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);

    // Creating mesh & supports:
    ddc::init_discrete_space<BSplinesR>(r_break_points);
    ddc::init_discrete_space<BSplinesTheta>(theta_break_points);

    NonUniformGridBase<R>::Impl<GridR, Kokkos::HostSpace> greville_sampling_r
            = SplineInterpPointsR::template get_sampling<GridR>();
    std::vector<CoordR> radial_points;
    Idx<GridR> ir = greville_sampling_r.front();
    for (int i(0); i < greville_sampling_r.size(); ++i, ++ir) {
        radial_points.push_back(greville_sampling_r.coordinate(ir));
    }
    radial_points[0] = CoordR(radial_points[1] * 0.5);

    ddc::init_discrete_space<GridR>(radial_points);
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
    const PhysicalToLogicalMapping to_logical_mapping;


    // --- Advection operator -------------------------------------------------------------------------
    SplineInterpolatorRTheta interpolator(grid);

    RK3Builder const time_stepper;
    PolarFootFinder find_feet = make_polar_foot_finder<TestFixture::FFSpace, TestFixture::AFSpace>(
            time_stepper,
            to_physical_mapping,
            grid,
            builder,
            spline_evaluator_extrapol);

    BslAdvectionPolar advection_operator(interpolator, find_feet, to_physical_mapping);

    // --- Advection field finder ---------------------------------------------------------------------
    AdvectionFieldFinder advection_field_computer(to_physical_mapping);


    // --- Choice of the simulation -------------------------------------------------------------------
    AdvectionFieldSimulation simulation
            = get_simulation<TestFixture::simulation_choice>(to_physical_mapping, r_min, r_max);

    // ================================================================================================
    // SIMULATION DATA                                                                                 |
    // ================================================================================================

    // ================================================================================================
    // INITIALISATION                                                                                 |
    // ================================================================================================
    host_t<DFieldMemRTheta> density_rtheta_alloc(grid);
    host_t<DFieldMemRTheta> density_xy_alloc(grid);

    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_exact_alloc(grid);
    host_t<DVectorFieldMemRTheta<R, Theta>> advection_field_rtheta_alloc(grid);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_alloc(grid);
    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_from_rtheta_alloc(grid);
    DVector<X, Y> advection_field_xy_centre;

    host_t<DFieldMemRTheta> electrostatic_potential_alloc(grid);

    host_t<DFieldRTheta> density_rtheta(density_rtheta_alloc);
    host_t<DFieldRTheta> density_xy(density_xy_alloc);
    host_t<DFieldRTheta> electrostatic_potential(electrostatic_potential_alloc);

    host_t<DVectorFieldRTheta<X, Y>> advection_field_exact(advection_field_exact_alloc);
    host_t<DVectorFieldRTheta<R, Theta>> advection_field_rtheta(advection_field_rtheta_alloc);
    host_t<DVectorFieldRTheta<X, Y>> advection_field_xy(advection_field_xy_alloc);
    host_t<DVectorFieldRTheta<X, Y>> advection_field_xy_from_rtheta(
            advection_field_xy_from_rtheta_alloc);



    // Initialise functions ******************************************
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));
        CoordXY const coord_xy(to_physical_mapping(coord_rtheta));

        density_rtheta(irtheta) = simulation.function(coord_rtheta);
        density_xy(irtheta) = density_rtheta(irtheta);
        electrostatic_potential(irtheta) = simulation.electrostatical_potential(coord_xy, 0);

        ddcHelper::assign_vector_field_element(
                advection_field_exact,
                irtheta,
                simulation.advection_field(coord_xy, 0));
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
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        ddcHelper::assign_vector_field_element(
                get_field(difference_between_fields_exact_and_xy),
                irtheta,
                advection_field_exact(irtheta) - advection_field_xy(irtheta));
    });


    // > Compare the advection field computed on RTheta to the advection field computed on XY
    host_t<DVectorFieldMemRTheta<X, Y>> difference_between_fields_xy_and_rtheta(grid);

    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

        // Jacobian matrix
        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> J
                = to_physical_mapping.jacobian_matrix(coord_rtheta);

        // computation made in BslAdvectionRTheta operator:
        ddcHelper::assign_vector_field_element(
                advection_field_xy_from_rtheta,
                irtheta,
                tensor_mul(index<'i', 'j'>(J), index<'j'>(advection_field_rtheta(irtheta))));

        // compare
        ddcHelper::assign_vector_field_element(
                get_field(difference_between_fields_xy_and_rtheta),
                irtheta,
                advection_field_xy_from_rtheta(irtheta) - advection_field_xy(irtheta));
    });

    // --- Check the difference on advection fields  --------------------------------------------------
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        EXPECT_LE(
                std::abs(ddcHelper::get<X>(difference_between_fields_exact_and_xy)(irtheta)),
                1e-5);
        EXPECT_LE(
                std::abs(ddcHelper::get<Y>(difference_between_fields_exact_and_xy)(irtheta)),
                1e-5);

        EXPECT_LE(
                std::abs(ddcHelper::get<X>(difference_between_fields_xy_and_rtheta)(irtheta)),
                1e-13);
        EXPECT_LE(
                std::abs(ddcHelper::get<Y>(difference_between_fields_xy_and_rtheta)(irtheta)),
                1e-13);
    });

    if constexpr (TestFixture::AFSpace == AdvectionFieldSpace::LOGICAL) {
        auto density_rtheta_device
                = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), density_rtheta);
        auto advection_field_rtheta_device = ddcHelper::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                advection_field_rtheta);

        // ================================================================================================
        // SIMULATION                                                                                     |
        // ================================================================================================
        for (int iter(0); iter < iter_nb; ++iter) {
            // --- operator() with logical advection field.
            advection_operator(
                    get_field(density_rtheta_device),
                    get_const_field(advection_field_rtheta_device),
                    dt);

            ddc::parallel_deepcopy(density_rtheta, get_const_field(density_rtheta_device));

            ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
                CoordXY foot = simulation.electrostatical_potential.exact_feet(
                        to_physical_mapping(ddc::coordinate(irtheta)),
                        dt * (iter + 1));
                double exact = simulation.function(to_logical_mapping(foot));
                EXPECT_NEAR(density_rtheta(irtheta), exact, 5e-7);
            });
        }

    } else {
        auto density_xy_device
                = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), density_xy);
        auto advection_field_xy_device = ddcHelper::
                create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), advection_field_xy);

        // ================================================================================================
        // SIMULATION                                                                                     |
        // ================================================================================================
        for (int iter(0); iter < iter_nb; ++iter) {
            // --- operator() with logical advection field.
            advection_operator(
                    get_field(density_xy_device),
                    get_const_field(advection_field_xy_device),
                    dt);

            ddc::parallel_deepcopy(density_xy, get_const_field(density_xy_device));

            ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
                CoordXY foot = simulation.electrostatical_potential.exact_feet(
                        to_physical_mapping(ddc::coordinate(irtheta)),
                        dt * (iter + 1));
                double exact = simulation.function(to_logical_mapping(foot));
                EXPECT_NEAR(density_xy(irtheta), exact, 5e-3);
            });
        }
    }

    end_simulation = std::chrono::system_clock::now();
    display_time_difference("Simulation time: ", start_simulation, end_simulation);
}
