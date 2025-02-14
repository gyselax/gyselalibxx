//#pragma once
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

#include "advection_simulation_utils.hpp"
#include "bsl_advection_rp.hpp"
#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "directional_tag.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "spline_polar_foot_finder.hpp"
#include "test_cases.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"


namespace {

namespace fs = std::filesystem;

using CircularToCartMapping = CircularToCartesian<R, Theta, X, Y>;
using CzarnyToCartMapping = CzarnyToCartesian<R, Theta, X, Y>;
using CartToCircularMapping = CartesianToCircular<X, Y, R, Theta>;
using CartToCzarnyMapping = CartesianToCzarny<X, Y, R, Theta>;
using CircularToPseudoCartMapping = CircularToCartesian<R, Theta, X_pC, Y_pC>;
using DiscreteMappingBuilder = DiscreteToCartesianBuilder<
        X,
        Y,
        SplineRThetaBuilder,
        SplineRThetaEvaluatorConstBound>;


} // end namespace
template <
        class LogicalToPhysicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class AnalyticalPhysicalToLogicalMapping,
        class AnalyticalLogicalToPhysicalMapping>
struct SimulationParameters
{
public:
    using X_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    using Y_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;

public:
    LogicalToPhysicalMapping const& to_physical_mapping;
    AnalyticalPhysicalToLogicalMapping const to_logical_mapping;
    LogicalToPseudoPhysicalMapping const& analytical_to_pseudo_physical_mapping;
    AnalyticalLogicalToPhysicalMapping const& analytical_to_physical_mapping;
    std::string mapping_name;
    std::string domain_name;
    SimulationParameters(
            LogicalToPhysicalMapping const& map,
            LogicalToPseudoPhysicalMapping const& pseudo_cart_map,
            AnalyticalPhysicalToLogicalMapping const& rev_map,
            AnalyticalLogicalToPhysicalMapping const& a_map,
            std::string m_name,
            std::string dom_name)
        : to_physical_mapping(map)
        , to_logical_mapping(rev_map)
        , analytical_to_pseudo_physical_mapping(pseudo_cart_map)
        , analytical_to_physical_mapping(a_map)
        , mapping_name(m_name)
        , domain_name(dom_name)
    {
    }
};

template <class TimeStepper>
struct NumericalMethodParameters
{
    TimeStepper time_stepper;
    double time_step;
    std::string method_name;
    NumericalMethodParameters(TimeStepper&& time_stepper, double step, std::string name)
        : time_stepper(std::move(time_stepper))
        , time_step(step)
        , method_name(name)
    {
    }
};


struct NumericalParams
{
    IdxRangeRTheta const grid;
    double const dt;

    NumericalParams(IdxRangeRTheta grid, double dt) : grid(grid), dt(dt) {};
    NumericalParams(NumericalParams&& params) = default;
    NumericalParams(NumericalParams& params) = default;
};


template <class X_adv, class Y_adv>
struct Numerics
{
private:
    NumericalParams params;

public:
    using ValFieldMem = FieldMemRTheta<CoordRTheta>;
    using DerivFieldMem = DVectorFieldMemRTheta<X_adv, Y_adv>;

    using NumericalTuple = std::tuple<
            NumericalMethodParameters<
                    Euler<ValFieldMem, DerivFieldMem, Kokkos::DefaultExecutionSpace>>,
            NumericalMethodParameters<
                    CrankNicolson<ValFieldMem, DerivFieldMem, Kokkos::DefaultExecutionSpace>>,
            NumericalMethodParameters<
                    RK3<ValFieldMem, DerivFieldMem, Kokkos::DefaultExecutionSpace>>,
            NumericalMethodParameters<
                    RK4<ValFieldMem, DerivFieldMem, Kokkos::DefaultExecutionSpace>>>;

    static constexpr int size_tuple = std::tuple_size<NumericalTuple> {};

    NumericalTuple numerics;

    explicit Numerics(NumericalParams m_params)
        : params(m_params)
        , numerics(std::make_tuple(
                  NumericalMethodParameters(
                          Euler<ValFieldMem, DerivFieldMem, Kokkos::DefaultExecutionSpace>(
                                  params.grid),
                          params.dt * 0.1,
                          "EULER"),
                  NumericalMethodParameters(
                          CrankNicolson<
                                  ValFieldMem,
                                  DerivFieldMem,
                                  Kokkos::DefaultExecutionSpace>(params.grid, 20, 1e-12),
                          params.dt,
                          "CRANK NICOLSON"),
                  NumericalMethodParameters(
                          RK3<ValFieldMem, DerivFieldMem, Kokkos::DefaultExecutionSpace>(
                                  params.grid),
                          params.dt,
                          "RK3"),
                  NumericalMethodParameters(
                          RK4<ValFieldMem, DerivFieldMem, Kokkos::DefaultExecutionSpace>(
                                  params.grid),
                          params.dt,
                          "RK4")))
    {
    }
};



struct GeneralParameters
{
    IdxRangeRTheta grid;
    PreallocatableSplineInterpolatorRTheta<ddc::NullExtrapolationRule> const& interpolator;
    SplineRThetaBuilder const& advection_builder;
    SplineRThetaEvaluatorConstBound& advection_evaluator;
    double final_time;
    bool if_save_curves;
    bool if_save_feet;
};

template <int i_map = 0, int i_feet = 0, class SimulationTuple>
void run_simulations_with_methods(
        SimulationTuple const& simulations,
        NumericalParams& num_params,
        GeneralParameters params)
{
    auto& sim = std::get<i_map>(simulations);

    using X_adv = typename std::remove_const_t<std::remove_reference_t<decltype(sim)>>::X_adv;
    using Y_adv = typename std::remove_const_t<std::remove_reference_t<decltype(sim)>>::Y_adv;

    Numerics<X_adv, Y_adv> methods(num_params);
    auto& num = std::get<i_feet>(methods.numerics);

    std::ostringstream name_stream;
    name_stream << sim.mapping_name << " MAPPING - " << sim.domain_name << " DOMAIN - "
                << num.method_name << " - ";
    std::string simulation_name = name_stream.str();

    std::ostringstream output_stream;
    output_stream << to_lower(sim.mapping_name) << "_" << to_lower(sim.domain_name) << "-"
                  << to_lower(num.method_name) << "-";
    std::string output_stem = output_stream.str();

    simulate_the_3_simulations(
            sim.to_physical_mapping,
            sim.to_logical_mapping,
            sim.analytical_to_pseudo_physical_mapping,
            sim.analytical_to_physical_mapping,
            params.grid,
            num.time_stepper,
            params.interpolator,
            params.advection_builder,
            params.advection_evaluator,
            params.final_time,
            num.time_step,
            params.if_save_curves,
            params.if_save_feet,
            output_stem,
            simulation_name);

    if constexpr (i_feet < methods.size_tuple - 1) {
        // Loop over numerical methods
        run_simulations_with_methods<i_map, i_feet + 1>(simulations, num_params, params);
    } else if constexpr (i_map < std::tuple_size<SimulationTuple> {} - 1) {
        // Loop over simulation parameters
        run_simulations_with_methods<i_map + 1, 0>(simulations, num_params, params);
    }
}

int main(int argc, char** argv)
{
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);


    PC_tree_t conf_voicexx;
    if (argc == 2) {
        conf_voicexx = PC_parse_path(fs::path(argv[1]).c_str());
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
    double const dt = PCpp_double(conf_voicexx, ".Time.time_step");
    double const final_time = PCpp_double(conf_voicexx, ".Time.final_time");
    bool const if_save_curves = PCpp_bool(conf_voicexx, ".Output.save_curves");
    bool const if_save_feet = PCpp_bool(conf_voicexx, ".Output.save_feet");

    // BUILD GRIDS ------------------------------------------------------------------------------
    // Grid creation of space. ------------------------------------------------------------------
    CoordTheta const p_min(0.0);
    CoordTheta const p_max(2.0 * M_PI);
    IdxStepTheta const p_ncells(PCpp_int(conf_voicexx, ".SplineMesh.p_ncells"));

    IdxRangeR const interpolation_idx_range_R = init_pseudo_uniform_spline_dependent_idx_range<
            GridR,
            BSplinesR,
            SplineInterpPointsR>(conf_voicexx, "r");
    PC_tree_destroy(&conf_voicexx);

    std::vector<CoordTheta> p_knots = build_uniform_break_points(p_min, p_max, p_ncells);
    ddc::init_discrete_space<BSplinesTheta>(p_knots);
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeTheta const interpolation_idx_range_P(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const grid(interpolation_idx_range_R, interpolation_idx_range_P);

    CoordR const rmin = ddc::coordinate(interpolation_idx_range_R.front());
    CoordR const rmax = ddc::coordinate(interpolation_idx_range_R.back());

    std::cout << "TESTS ON THE ADVECTION OPERATOR "
              << "FOR [rmin, rmax] = [" << double(rmin) << ", " << double(rmax) << "], "
              << "WITH NrxNt = " << interpolation_idx_range_R.size() << "x"
              << interpolation_idx_range_P.size() << " AND dt = " << dt << ": " << std::endl;



    // DEFINITION OF OPERATORS ------------------------------------------------------------------
    // --- Builders for the test function and the to_physical_mapping:
    SplineRThetaBuilder const builder(grid);
    SplineRThetaBuilder_host const builder_host(grid);

    // --- Evaluator for the test function:
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> p_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);
    SplineRThetaEvaluatorNullBound_host spline_evaluator_host(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    PreallocatableSplineInterpolatorRTheta interpolator(builder_host, spline_evaluator_host);


    // --- Evaluator for the test advection field:
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(rmin);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(rmax);


    SplineRThetaEvaluatorConstBound spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());



    // SET THE DIFFERENT PARAMETERS OF THE TESTS ------------------------------------------------
    CircularToCartMapping const from_circ_map;
    CircularToPseudoCartMapping const to_pseudo_circ_map;
    CartesianToCircular<X, Y, R, Theta> to_circ_map;
    CzarnyToCartMapping const from_czarny_map(0.3, 1.4);
    CartesianToCzarny<X, Y, R, Theta> const to_czarny_map(0.3, 1.4);
    DiscreteMappingBuilder const discrete_czarny_map_builder(
            Kokkos::DefaultExecutionSpace(),
            from_czarny_map,
            builder,
            spline_evaluator_extrapol);
    DiscreteToCartesian const from_discrete_czarny_map = discrete_czarny_map_builder();

    std::tuple simulations = std::make_tuple(
            SimulationParameters(
                    from_circ_map,
                    from_circ_map,
                    to_circ_map,
                    from_circ_map,
                    "CIRCULAR",
                    "PHYSICAL"),
            SimulationParameters(
                    from_czarny_map,
                    from_czarny_map,
                    to_czarny_map,
                    from_czarny_map,
                    "CZARNY",
                    "PHYSICAL"),
            SimulationParameters(
                    from_czarny_map,
                    to_pseudo_circ_map,
                    to_czarny_map,
                    from_czarny_map,
                    "CZARNY",
                    "PSEUDO CARTESIAN"),
            SimulationParameters(
                    from_discrete_czarny_map,
                    to_pseudo_circ_map,
                    to_czarny_map,
                    from_czarny_map,
                    "DISCRETE",
                    "PSEUDO CARTESIAN"));


    NumericalParams num_params = {grid, dt};


    // TO CLOCK THE SIMULATION --------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> start_full_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_full_simulation;

    start_full_simulation = std::chrono::system_clock::now();


    // SIMULATION: ==========================================================================

    GeneralParameters params
            = {grid,
               interpolator,
               builder,
               spline_evaluator_extrapol,
               final_time,
               if_save_curves,
               if_save_feet};
    run_simulations_with_methods(simulations, num_params, params);

    end_full_simulation = std::chrono::system_clock::now();
    display_time(start_full_simulation, end_full_simulation, "   Full simulation time:    ");
}
