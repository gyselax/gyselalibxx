#include <array>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <tuple>
#include <typeinfo>
#include <vector>

#include <ddc/ddc.hpp>

#include "../../advection/r_theta_test_cases.hpp"

#include "advection_simulation_utils.hpp"
#include "bsl_advection_polar.hpp"
#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "discrete_poloidal_cs_spline_mapping.hpp"
#include "discrete_poloidal_cs_spline_mapping_builder.hpp"
#include "euler.hpp"
#include "geometry_r_theta.hpp"
#include "input.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polar_foot_finder.hpp"
#include "polar_spline_evaluator.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "spline_definitions_r_theta.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"


namespace {

namespace fs = std::filesystem;

using CircularToCartMapping = CircularToCartesian<R, Theta, X, Y>;
using CzarnyToCartMapping = CzarnyToCartesian<R, Theta, X, Y>;
using CartToCircularMapping = CartesianToCircular<X, Y, R, Theta>;
using CartToCzarnyMapping = CartesianToCzarny<X, Y, R, Theta>;
using CircularToPseudoCartMapping = CircularToCartesian<R, Theta, X_pC, Y_pC>;
using DiscreteMappingBuilderHost = DiscretePoloidalCSSplineMappingBuilder<
        X,
        Y,
        SplineRThetaBuilder_host,
        SplineRThetaEvaluatorConstBound_host>;
using DiscreteMappingBuilder = DiscretePoloidalCSSplineMappingBuilder<
        X,
        Y,
        SplineRThetaBuilder,
        SplineRThetaEvaluatorConstBound>;


} // end namespace

enum TimeStepperChoice { TS_EULER, TS_CRANK_NICOLSON, TS_RK3, TS_RK4 };

template <TimeStepperChoice TSChoice>
auto get_time_stepper_builder()
{
    if constexpr (TSChoice == TS_EULER) {
        return EulerBuilder();
    } else if constexpr (TSChoice == TS_CRANK_NICOLSON) {
        return CrankNicolsonBuilder(20, 1e-12);
    } else if constexpr (TSChoice == TS_RK3) {
        return RK3Builder();
    } else if constexpr (TSChoice == TS_RK4) {
        return RK4Builder();
    }
}

struct GeneralParameters
{
    IdxRangeRTheta grid;
    SplineInterpolatorRTheta const& interpolator;
    SplineRThetaBuilder const& advection_builder;
    SplineRThetaEvaluatorConstBound& advection_evaluator;
    double final_time;
    bool if_save_curves;
    bool if_save_feet;
};

template <
        TimeStepperChoice TSChoice,
        class LogicalToPhysicalMappingHost,
        class LogicalToPhysicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class AnalyticalPhysicalToLogicalMapping,
        class AnalyticalLogicalToPhysicalMapping>
void run_simulations_with_foot_finder_method(
        double const dt,
        GeneralParameters params,
        LogicalToPhysicalMappingHost const& to_physical_mapping_host,
        LogicalToPhysicalMapping const& to_physical_mapping,
        LogicalToPseudoPhysicalMapping const& analytical_to_pseudo_physical_mapping,
        AnalyticalPhysicalToLogicalMapping const& to_logical_mapping,
        AnalyticalLogicalToPhysicalMapping const& analytical_to_physical_mapping,
        std::string const& simulation_name,
        std::string const& output_stem)
{
    constexpr FootFindingSpace FFSpace
            = std::is_same_v<
                      typename LogicalToPseudoPhysicalMapping::CoordResult,
                      Coord<X_pC, Y_pC>>
                      ? FootFindingSpace::PSEUDO_PHYSICAL
                      : FootFindingSpace::PHYSICAL;
    PolarFootFinder foot_finder = make_polar_foot_finder<FFSpace, AdvectionFieldSpace::PHYSICAL>(
            get_time_stepper_builder<TSChoice>(),
            to_physical_mapping,
            params.grid,
            params.advection_builder,
            params.advection_evaluator);

    BslAdvectionPolar advection_operator(params.interpolator, foot_finder, to_physical_mapping);

    run_simulations(
            to_physical_mapping_host,
            to_physical_mapping,
            to_logical_mapping,
            analytical_to_physical_mapping,
            params.grid,
            foot_finder,
            advection_operator,
            params.final_time,
            dt,
            params.if_save_curves,
            params.if_save_feet,
            output_stem,
            simulation_name);
}

template <
        class LogicalToPhysicalMappingHost,
        class LogicalToPhysicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class AnalyticalPhysicalToLogicalMapping,
        class AnalyticalLogicalToPhysicalMapping>
void run_simulations_with_methods(
        double const dt,
        GeneralParameters params,
        LogicalToPhysicalMappingHost const& to_physical_mapping_host,
        LogicalToPhysicalMapping const& to_physical_mapping,
        LogicalToPseudoPhysicalMapping const& analytical_to_pseudo_physical_mapping,
        AnalyticalPhysicalToLogicalMapping const& to_logical_mapping,
        AnalyticalLogicalToPhysicalMapping const& analytical_to_physical_mapping,
        std::string const& mapping_name,
        std::string const& domain_name)
{
    std::string simulation_name = mapping_name + " MAPPING - " + domain_name + " DOMAIN -";
    std::string output_stem = to_lower(mapping_name) + "_" + to_lower(domain_name) + "-";

    run_simulations_with_foot_finder_method<TS_EULER>(
            dt * 0.1,
            params,
            to_physical_mapping_host,
            to_physical_mapping,
            analytical_to_pseudo_physical_mapping,
            to_logical_mapping,
            analytical_to_physical_mapping,
            simulation_name + "EULER - ",
            output_stem + "euler-");

    run_simulations_with_foot_finder_method<TS_CRANK_NICOLSON>(
            dt,
            params,
            to_physical_mapping_host,
            to_physical_mapping,
            analytical_to_pseudo_physical_mapping,
            to_logical_mapping,
            analytical_to_physical_mapping,
            simulation_name + "CRANK NICOLSON - ",
            output_stem + "crank nicolson-");

    run_simulations_with_foot_finder_method<TS_RK3>(
            dt,
            params,
            to_physical_mapping_host,
            to_physical_mapping,
            analytical_to_pseudo_physical_mapping,
            to_logical_mapping,
            analytical_to_physical_mapping,
            simulation_name + "RK3 - ",
            output_stem + "rk3-");

    run_simulations_with_foot_finder_method<TS_RK4>(
            dt,
            params,
            to_physical_mapping_host,
            to_physical_mapping,
            analytical_to_pseudo_physical_mapping,
            to_logical_mapping,
            analytical_to_physical_mapping,
            simulation_name + "RK4 - ",
            output_stem + "rk4-");
}

int main(int argc, char** argv)
{
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
    PC_errhandler(PC_NULL_HANDLER);



    // Parameters of the simulation. ------------------------------------------------------------
    double const dt = PCpp_double(conf_gyselalibxx, ".Time.time_step");
    double const final_time = PCpp_double(conf_gyselalibxx, ".Time.final_time");
    bool const if_save_curves = PCpp_bool(conf_gyselalibxx, ".Output.save_curves");
    bool const if_save_feet = PCpp_bool(conf_gyselalibxx, ".Output.save_feet");

    // BUILD GRIDS ------------------------------------------------------------------------------
    // Grid creation of space. ------------------------------------------------------------------
    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(PCpp_int(conf_gyselalibxx, ".SplineMesh.theta_ncells"));

    IdxRangeR const interpolation_idx_range_r = init_pseudo_uniform_spline_dependent_idx_range<
            GridR,
            BSplinesR,
            SplineInterpPointsR>(conf_gyselalibxx, "r");
    PC_tree_destroy(&conf_gyselalibxx);

    std::vector<CoordTheta> theta_break_points
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);
    ddc::init_discrete_space<BSplinesTheta>(theta_break_points);
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeTheta const interpolation_idx_range_theta(
            SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const grid(interpolation_idx_range_r, interpolation_idx_range_theta);

    CoordR const rmin = ddc::coordinate(interpolation_idx_range_r.front());
    CoordR const rmax = ddc::coordinate(interpolation_idx_range_r.back());

    std::cout << "TESTS ON THE ADVECTION OPERATOR "
              << "FOR [rmin, rmax] = [" << double(rmin) << ", " << double(rmax) << "], "
              << "WITH NrxNt = " << interpolation_idx_range_r.size() << "x"
              << interpolation_idx_range_theta.size() << " AND dt = " << dt << ": " << std::endl;



    // DEFINITION OF OPERATORS ------------------------------------------------------------------
    // --- Builders for the test function and the to_physical_mapping:
    SplineRThetaBuilder_host const builder_host(grid);
    SplineRThetaBuilder const builder(grid);

    // --- Interpolator for the test function:
    SplineInterpolatorRTheta spline_interpolator(grid);

    // --- Evaluator for the test advection field:
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(rmin);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(rmax);


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



    // SET THE DIFFERENT PARAMETERS OF THE TESTS ------------------------------------------------
    // Offset the centre of the circle to ensure that this is correctly handled
    Coord<X, Y> origin_point(6.2, 0.8);
    CircularToCartMapping const from_circ_map(origin_point);
    Coord<X_pC, Y_pC> origin_point_pc(6.2, 0.8);
    CircularToPseudoCartMapping const to_pseudo_circ_map(origin_point_pc);
    CartesianToCircular<X, Y, R, Theta> to_circ_map(origin_point);
    CzarnyToCartMapping const from_czarny_map(0.3, 1.4, origin_point);
    CartesianToCzarny<X, Y, R, Theta> const to_czarny_map(0.3, 1.4, origin_point);
    DiscreteMappingBuilderHost const discrete_czarny_map_builder_host(
            Kokkos::DefaultHostExecutionSpace(),
            from_czarny_map,
            builder_host,
            spline_evaluator_extrapol_host);
    DiscretePoloidalCSSplineMapping const from_discrete_czarny_map_host
            = discrete_czarny_map_builder_host();
    DiscreteMappingBuilder const discrete_czarny_map_builder(
            Kokkos::DefaultExecutionSpace(),
            from_czarny_map,
            builder,
            spline_evaluator_extrapol);
    DiscretePoloidalCSSplineMapping const from_discrete_czarny_map = discrete_czarny_map_builder();


    // TO CLOCK THE SIMULATION --------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> start_full_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_full_simulation;

    start_full_simulation = std::chrono::system_clock::now();


    // SIMULATION: ==========================================================================

    GeneralParameters params
            = {grid,
               spline_interpolator,
               builder,
               spline_evaluator_extrapol,
               final_time,
               if_save_curves,
               if_save_feet};

    run_simulations_with_methods(
            dt,
            params,
            from_circ_map,
            from_circ_map,
            from_circ_map,
            to_circ_map,
            from_circ_map,
            "CIRCULAR",
            "PHYSICAL");
    run_simulations_with_methods(
            dt,
            params,
            from_czarny_map,
            from_czarny_map,
            from_czarny_map,
            to_czarny_map,
            from_czarny_map,
            "CZARNY",
            "PHYSICAL");
    run_simulations_with_methods(
            dt,
            params,
            from_czarny_map,
            from_czarny_map,
            to_pseudo_circ_map,
            to_czarny_map,
            from_czarny_map,
            "CZARNY",
            "PSEUDO CARTESIAN");
    run_simulations_with_methods(
            dt,
            params,
            from_discrete_czarny_map_host,
            from_discrete_czarny_map,
            to_pseudo_circ_map,
            to_czarny_map,
            from_czarny_map,
            "DISCRETE",
            "PSEUDO CARTESIAN");

    end_full_simulation = std::chrono::system_clock::now();
    display_time(start_full_simulation, end_full_simulation, "   Full simulation time:    ");
}
