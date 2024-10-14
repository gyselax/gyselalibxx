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

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_builder.hpp>
#include <sll/mapping/discrete_to_cartesian.hpp>
#include <sll/math_tools.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>

#include "advection_domain.hpp"
#include "advection_simulation_utils.hpp"
#include "bsl_advection_rp.hpp"
#include "crank_nicolson.hpp"
#include "ddc_aliases.hpp"
#include "directional_tag.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "mesh_builder.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "test_cases.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"


namespace {

namespace fs = std::filesystem;

using CurvilinearMapping = Curvilinear2DToCartesian<X, Y, R, Theta>;
using CircularMapping = CircularToCartesian<X, Y, R, Theta>;
using CzarnyMapping = CzarnyToCartesian<X, Y, R, Theta>;
using DiscreteMappingBuilder
        = DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>;


} // end namespace
template <class Mapping, class AnalyticalMapping, class AdvectionDomain>
struct SimulationParameters
{
public:
    Mapping const& mapping;
    AnalyticalMapping const& analytical_mapping;
    AdvectionDomain const& advection_domain;
    std::string mapping_name;
    std::string domain_name;
    using IdxRangeSimulationAdvection = AdvectionDomain;
    SimulationParameters(
            Mapping const& map,
            AnalyticalMapping const& a_map,
            AdvectionDomain const& advection_dom,
            std::string m_name,
            std::string dom_name)
        : mapping(map)
        , analytical_mapping(a_map)
        , advection_domain(advection_dom)
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
    IdxRangeRTheta grid;
    double dt;

    NumericalParams(IdxRangeRTheta grid, double dt) : grid(grid), dt(dt) {};
    NumericalParams(NumericalParams&& params) = default;
    NumericalParams(NumericalParams& params) = default;
};


template <class AdvectionDomain>
struct Numerics
{
private:
    AdvectionDomain advection_domain;
    NumericalParams params;

public:
    using X_adv = typename AdvectionDomain::X_adv;
    using Y_adv = typename AdvectionDomain::Y_adv;

    using ValFieldMem = host_t<FieldMemRTheta<CoordRTheta>>;
    using DerivFieldMem = host_t<DVectorFieldMemRTheta<X_adv, Y_adv>>;

    using NumericalTuple = std::tuple<
            NumericalMethodParameters<Euler<ValFieldMem, DerivFieldMem>>,
            NumericalMethodParameters<CrankNicolson<ValFieldMem, DerivFieldMem>>,
            NumericalMethodParameters<RK3<ValFieldMem, DerivFieldMem>>,
            NumericalMethodParameters<RK4<ValFieldMem, DerivFieldMem>>>;

    static constexpr int size_tuple = std::tuple_size<NumericalTuple> {};

    NumericalTuple numerics;

    Numerics(AdvectionDomain m_advection_domain, NumericalParams m_params)
        : advection_domain(m_advection_domain)
        , params(m_params)
        , numerics(std::make_tuple(
                  NumericalMethodParameters(
                          Euler<ValFieldMem, DerivFieldMem>(params.grid),
                          params.dt * 0.1,
                          "EULER"),
                  NumericalMethodParameters(
                          CrankNicolson<ValFieldMem, DerivFieldMem>(params.grid, 20, 1e-12),
                          params.dt,
                          "CRANK NICOLSON"),
                  NumericalMethodParameters(
                          RK3<ValFieldMem, DerivFieldMem>(params.grid),
                          params.dt,
                          "RK3"),
                  NumericalMethodParameters(
                          RK4<ValFieldMem, DerivFieldMem>(params.grid),
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

    Numerics methods(sim.advection_domain, num_params);
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
            sim.mapping,
            sim.analytical_mapping,
            params.grid,
            num.time_stepper,
            sim.advection_domain,
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
    // --- Builders for the test function and the mapping:
    SplineRThetaBuilder const builder(grid);

    // --- Evaluator for the test function:
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> p_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    PreallocatableSplineInterpolatorRTheta interpolator(builder, spline_evaluator);


    // --- Evaluator for the test advection field:
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(rmin);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(rmax);


    SplineRThetaEvaluatorConstBound spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());



    // SET THE DIFFERENT PARAMETERS OF THE TESTS ------------------------------------------------
    CircularMapping const circ_map;
    CzarnyMapping const czarny_map(0.3, 1.4);
    DiscreteMappingBuilder const discrete_czarny_map_builder(
            Kokkos::DefaultHostExecutionSpace(),
            czarny_map,
            builder,
            spline_evaluator_extrapol);
    DiscreteToCartesian const discrete_czarny_map = discrete_czarny_map_builder();

    AdvectionPhysicalDomain<CircularMapping> const physical_circular_mapping(circ_map);
    AdvectionPhysicalDomain<CzarnyMapping> const physical_czarny_mapping(czarny_map);
    AdvectionPseudoCartesianDomain<CzarnyMapping> const pseudo_cartesian_czarny_mapping(czarny_map);

    std::tuple simulations = std::make_tuple(
            SimulationParameters(
                    circ_map,
                    circ_map,
                    physical_circular_mapping,
                    "CIRCULAR",
                    "PHYSICAL"),
            SimulationParameters(
                    czarny_map,
                    czarny_map,
                    physical_czarny_mapping,
                    "CZARNY",
                    "PHYSICAL"),
            SimulationParameters(
                    czarny_map,
                    czarny_map,
                    pseudo_cartesian_czarny_mapping,
                    "CZARNY",
                    "PSEUDO CARTESIAN"),
            SimulationParameters(
                    discrete_czarny_map,
                    czarny_map,
                    pseudo_cartesian_czarny_mapping,
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
