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

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>
#include <sll/math_tools.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include <crank_nicolson.hpp>
#include <directional_tag.hpp>
#include <euler.hpp>
#include <rk2.hpp>
#include <rk3.hpp>
#include <rk4.hpp>
#include <species_info.hpp>
#include <stdio.h>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

#include "advection_domain.hpp"
#include "advection_simulation_utils.hpp"
#include "bsl_advection_rp.hpp"
#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "test_cases.hpp"


namespace {

namespace fs = std::filesystem;

using CurvilinearMapping = Curvilinear2DToCartesian<RDimX, RDimY, RDimR, RDimP>;
using AnalyticalInvertibleMapping
        = AnalyticalInvertibleCurvilinear2DToCartesian<RDimX, RDimY, RDimR, RDimP>;
using CircularMapping = CircularToCartesian<RDimX, RDimY, RDimR, RDimP>;
using CzarnyMapping = CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP>;
using DiscreteMapping = DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>;


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
    using SimultationAdvectionDomain = AdvectionDomain;
    SimulationParameters(
            Mapping const& map,
            AnalyticalMapping const& a_map,
            AdvectionDomain const& dom,
            std::string m_name,
            std::string dom_name)
        : mapping(map)
        , analytical_mapping(a_map)
        , advection_domain(dom)
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
    IDomainRP grid;
    double dt;

    NumericalParams(IDomainRP grid, double dt) : grid(grid), dt(dt) {};
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
    using RDimX_adv = typename AdvectionDomain::RDimX_adv;
    using RDimY_adv = typename AdvectionDomain::RDimY_adv;

    using ValChunk = FieldRP<CoordRP>;
    using DerivChunk = VectorDFieldRP<RDimX_adv, RDimY_adv>;

    using NumericalTuple = std::tuple<
            NumericalMethodParameters<Euler<ValChunk, DerivChunk>>,
            NumericalMethodParameters<CrankNicolson<ValChunk, DerivChunk>>,
            NumericalMethodParameters<RK3<ValChunk, DerivChunk>>,
            NumericalMethodParameters<RK4<ValChunk, DerivChunk>>>;

    static constexpr int size_tuple = std::tuple_size<NumericalTuple> {};

    NumericalTuple numerics;

    Numerics(AdvectionDomain m_advection_domain, NumericalParams m_params)
        : advection_domain(m_advection_domain)
        , params(m_params)
        , numerics(std::make_tuple(
                  NumericalMethodParameters(
                          Euler<ValChunk, DerivChunk>(params.grid),
                          params.dt * 0.1,
                          "EULER"),
                  NumericalMethodParameters(
                          CrankNicolson<ValChunk, DerivChunk>(params.grid, 20, 1e-12),
                          params.dt,
                          "CRANK NICOLSON"),
                  NumericalMethodParameters(
                          RK3<ValChunk, DerivChunk>(params.grid),
                          params.dt,
                          "RK3"),
                  NumericalMethodParameters(
                          RK4<ValChunk, DerivChunk>(params.grid),
                          params.dt,
                          "RK4")))
    {
    }
};



struct GeneralParameters
{
    IDomainRP grid;
    PreallocatableSplineInterpolatorRP const& interpolator;
    SplineRPBuilder const& advection_builder;
    SplineRPEvaluator& advection_evaluator;
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
    ::ddc::ScopeGuard scope(argc, argv);


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



    // Parameters of the grid. ---------------------------------------------------------------
    double const rmin = PCpp_double(conf_voicexx, ".Mesh.r_min");
    double const rmax = PCpp_double(conf_voicexx, ".Mesh.r_max");
    int const Nr = PCpp_int(conf_voicexx, ".Mesh.r_size");
    int const Nt = PCpp_int(conf_voicexx, ".Mesh.p_size");
    double const dt = PCpp_double(conf_voicexx, ".Mesh.time_step");
    double const final_time = PCpp_double(conf_voicexx, ".Mesh.final_time");
    bool const if_save_curves = PCpp_bool(conf_voicexx, ".Mesh.save_curves");
    bool const if_save_feet = PCpp_bool(conf_voicexx, ".Mesh.save_feet");
    PC_tree_destroy(&conf_voicexx);

    std::cout << "TESTS ON THE ADVECTION OPERATOR "
              << "FOR [rmin, rmax] = [" << rmin << ", " << rmax << "], "
              << "WITH NrxNt = " << Nr << "x" << Nt << " AND dt = " << dt << ": " << std::endl;

    // BUILD GRIDS ------------------------------------------------------------------------------
    // Grid creation of space. ------------------------------------------------------------------
    CoordR const r_min(rmin);
    CoordR const r_max(rmax);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordP> p_knots(p_size + 1);

    for (int i(0); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    r_knots[r_size] = CoordR(r_max);
    for (int i(0); i < p_size + 1; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<IDimR>(SplineInterpPointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(SplineInterpPointsP::get_sampling());

    IDomainR const interpolation_domain_R(SplineInterpPointsR::get_domain());
    IDomainP const interpolation_domain_P(SplineInterpPointsP::get_domain());
    IDomainRP const grid(interpolation_domain_R, interpolation_domain_P);



    // DEFINITION OF OPERATORS ------------------------------------------------------------------
    // --- Builders for the test function and the mapping:
    SplineRPBuilder const builder(grid);

    // --- Evaluator for the test function:
    SplineRPEvaluator spline_evaluator(
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>);

    PreallocatableSplineInterpolatorRP interpolator(builder, spline_evaluator);


    // --- Evaluator for the test advection field:
    ConstantExtrapolationBoundaryValue2D<BSplinesR, BSplinesP, RDimR> boundary_condition_r_left(
            r_min);
    ConstantExtrapolationBoundaryValue2D<BSplinesR, BSplinesP, RDimR> boundary_condition_r_right(
            r_max);


    SplineRPEvaluator spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>);



    // SET THE DIFFERENT PARAMETERS OF THE TESTS ------------------------------------------------
    CircularMapping const circ_map;
    CzarnyMapping const czarny_map(0.3, 1.4);
    DiscreteMapping const discrete_czarny_map(
            DiscreteMapping::
                    analytical_to_discrete(czarny_map, builder, spline_evaluator_extrapol));

    AdvectionPhysicalDomain<AnalyticalInvertibleMapping> const physical_circular_mapping(circ_map);
    AdvectionPhysicalDomain<AnalyticalInvertibleMapping> const physical_czarny_mapping(czarny_map);
    AdvectionPseudoCartesianDomain<CzarnyMapping> const pseudo_cartesian_czarny_mapping(czarny_map);

    auto simulations = std::make_tuple(
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
