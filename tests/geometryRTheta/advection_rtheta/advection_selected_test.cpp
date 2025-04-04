
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <typeinfo>

#include <ddc/ddc.hpp>

#include "advection_simulation_utils.hpp"
#include "bsl_advection_polar.hpp"
#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_helper.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "input.hpp"
#include "itimestepper.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "spline_polar_foot_finder.hpp"
#include "test_cases.hpp"



namespace fs = std::filesystem;

namespace {
#if defined(CIRCULAR_MAPPING_PHYSICAL)
using X_adv = X;
using Y_adv = Y;
#elif defined(CZARNY_MAPPING_PHYSICAL)
using X_adv = X;
using Y_adv = Y;

#elif defined(CZARNY_MAPPING_PSEUDO_CARTESIAN)
using X_adv = X_pC;
using Y_adv = Y_pC;

#elif defined(DISCRETE_MAPPING_PSEUDO_CARTESIAN)
using X_adv = X_pC;
using Y_adv = Y_pC;
#endif

} //end namespace

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



    // Parameters of the grid. ---------------------------------------------------------------
    double const dt = PCpp_double(conf_gyselalibxx, ".Time.time_step");
    double const final_time = PCpp_double(conf_gyselalibxx, ".Time.final_time");
    bool const save_curves = PCpp_bool(conf_gyselalibxx, ".Output.save_curves");
    bool const save_feet = PCpp_bool(conf_gyselalibxx, ".Output.save_feet");

    if (save_curves or save_feet) {
        fs::create_directory("output");
    }
    if (save_curves) {
        fs::create_directory("output/curves");
    }

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

    std::vector<CoordTheta> theta_knots
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);
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

    std::ofstream file("r_knots.txt");
    for_each(interpolation_idx_range_r, [&](IdxR ir) {
        file << (ir - interpolation_idx_range_r.front()).value() << " "
             << double(ddc::coordinate(ir)) << std::endl;
    });
    file.close();



    // DEFINITION OF OPERATORS ------------------------------------------------------------------
    // --- Builders for the test function and the to_physical_mapping:
    SplineRThetaBuilder_host const builder_host(grid);
    SplineRThetaBuilder const builder(grid);

    // --- Evaluator for the test function:
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

    PreallocatableSplineInterpolator2D interpolator(builder, spline_evaluator);


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


    std::string key;

    // SELECTION OF THE MAPPING AND THE ADVECTION DOMAIN.
#if defined(CIRCULAR_MAPPING_PHYSICAL)
    CircularToCartesian<R, Theta, X, Y> to_physical_analytical_mapping;
    CircularToCartesian<R, Theta, X, Y> to_physical_mapping;
    CircularToCartesian<R, Theta, X, Y> to_physical_mapping_host;
    CartesianToCircular<X, Y, R, Theta> to_logical_analytical_mapping;
    CircularToCartesian<R, Theta, X, Y> const& logical_to_pseudo_cart_mapping(
            to_physical_analytical_mapping);
    std::string const mapping_name = "CIRCULAR";
    std::string const adv_domain_name = "PHYSICAL";
    key += "circular_physical";
#else

    double const czarny_e = 0.3;
    double const czarny_epsilon = 1.4;

#if defined(CZARNY_MAPPING_PHYSICAL)
    CzarnyToCartesian<R, Theta, X, Y> to_physical_analytical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<R, Theta, X, Y> to_physical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<R, Theta, X, Y> to_physical_mapping_host(czarny_e, czarny_epsilon);
    CartesianToCzarny<X, Y, R, Theta> to_logical_analytical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<R, Theta, X, Y> const& logical_to_pseudo_cart_mapping(
            to_physical_analytical_mapping);
    std::string const mapping_name = "CZARNY";
    std::string const adv_domain_name = "PHYSICAL";
    key += "czarny_physical";

#elif defined(CZARNY_MAPPING_PSEUDO_CARTESIAN)
    CzarnyToCartesian<R, Theta, X, Y> to_physical_analytical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<R, Theta, X, Y> to_physical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<R, Theta, X, Y> to_physical_mapping_host(czarny_e, czarny_epsilon);
    CartesianToCzarny<X, Y, R, Theta> to_logical_analytical_mapping(czarny_e, czarny_epsilon);
    CircularToCartesian<R, Theta, X_pC, Y_pC> logical_to_pseudo_cart_mapping;
    std::string const mapping_name = "CZARNY";
    std::string const adv_domain_name = "PSEUDO CARTESIAN";
    key += "czarny_pseudo_cartesian";

#elif defined(DISCRETE_MAPPING_PSEUDO_CARTESIAN)
    CzarnyToCartesian<R, Theta, X, Y> to_physical_analytical_mapping(czarny_e, czarny_epsilon);
    CartesianToCzarny<X, Y, R, Theta> to_logical_analytical_mapping(czarny_e, czarny_epsilon);
    DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluatorConstBound_host>
            mapping_builder_host(
                    Kokkos::DefaultHostExecutionSpace(),
                    to_physical_analytical_mapping,
                    builder_host,
                    spline_evaluator_extrapol_host);
    DiscreteToCartesian to_physical_mapping_host = mapping_builder_host();
    DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>
            mapping_builder(
                    Kokkos::DefaultExecutionSpace(),
                    to_physical_analytical_mapping,
                    builder,
                    spline_evaluator_extrapol);
    DiscreteToCartesian to_physical_mapping = mapping_builder();
    CircularToCartesian<R, Theta, X_pC, Y_pC> logical_to_pseudo_cart_mapping;
    std::string const mapping_name = "DISCRETE";
    std::string const adv_domain_name = "PSEUDO CARTESIAN";
    key += "discrete_pseudo_cartesian";
#endif
#endif

    key += "-";

    // SELECTION OF THE TIME INTEGRATION METHOD.
#if defined(EULER_METHOD)
    Euler<FieldMemRTheta<CoordRTheta>,
          DVectorFieldMemRTheta<X_adv, Y_adv>,
          Kokkos::DefaultExecutionSpace>
            time_stepper(grid);
    std::string const method_name = "EULER";
    key += "euler";

#elif defined(CRANK_NICOLSON_METHOD)
    double const epsilon_CN = 1e-8;
    CrankNicolson<
            FieldMemRTheta<CoordRTheta>,
            DVectorFieldMemRTheta<X_adv, Y_adv>,
            Kokkos::DefaultExecutionSpace>
            time_stepper(grid, 20, epsilon_CN);
    std::string const method_name = "CRANK NICOLSON";
    key += "crank_nicolson";

#elif defined(RK3_METHOD)
    RK3<FieldMemRTheta<CoordRTheta>,
        DVectorFieldMemRTheta<X_adv, Y_adv>,
        Kokkos::DefaultExecutionSpace>
            time_stepper(grid);
    std::string const method_name = "RK3";
    key += "rk3";

#elif defined(RK4_METHOD)
    RK4<FieldMemRTheta<CoordRTheta>,
        DVectorFieldMemRTheta<X_adv, Y_adv>,
        Kokkos::DefaultExecutionSpace>
            time_stepper(grid);
    std::string const method_name = "RK4";
    key += "rk4";
#endif

    key += "-";

    // SELECTION OF THE SIMULATION.
#if defined(TRANSLATION_SIMULATION)
    AdvectionSimulation simulation
            = get_translation_simulation(to_physical_mapping_host, rmin, rmax);
    std::string const simu_type = "TRANSLATION";
    key += "Translation";

#elif defined(ROTATION_SIMULATION)
    AdvectionSimulation simulation = get_rotation_simulation(to_physical_mapping_host, rmin, rmax);
    std::string const simu_type = "ROTATION";
    key += "Rotation";

#elif defined(DECENTRED_ROTATION_SIMULATION)
    AdvectionSimulation simulation = get_decentred_rotation_simulation(to_physical_mapping_host);
    std::string const simu_type = "DECENTRED ROTATION";
    key += "Decentred_rotation";
#endif

    std::string output_folder = key + "_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }

    SplinePolarFootFinder const foot_finder(
            time_stepper,
            to_physical_mapping,
            logical_to_pseudo_cart_mapping,
            builder,
            spline_evaluator_extrapol);

    BslAdvectionPolar advection_operator(interpolator, foot_finder, to_physical_mapping);

    std::cout << mapping_name << " MAPPING - " << adv_domain_name << " DOMAIN - " << method_name
              << " - " << simu_type << " : " << std::endl;
    simulate(
            to_physical_mapping_host,
            to_physical_mapping,
            to_logical_analytical_mapping,
            logical_to_pseudo_cart_mapping,
            to_physical_analytical_mapping,
            grid,
            foot_finder,
            advection_operator,
            simulation,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);
}
