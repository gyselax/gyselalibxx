
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <typeinfo>

#include <ddc/ddc.hpp>

#include <sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>
#include <sll/math_tools.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>

#include "advection_domain.hpp"
#include "advection_simulation_utils.hpp"
#include "bsl_advection_rp.hpp"
#include "crank_nicolson.hpp"
#include "ddc_helper.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "itimestepper.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "spline_foot_finder.hpp"
#include "test_cases.hpp"



namespace fs = std::filesystem;

namespace {
#if defined(CIRCULAR_MAPPING_PHYSICAL)
using X_adv = typename AdvectionPhysicalDomain<CircularToCartesian<X, Y, R, Theta>>::X_adv;
using Y_adv = typename AdvectionPhysicalDomain<CircularToCartesian<X, Y, R, Theta>>::Y_adv;
#elif defined(CZARNY_MAPPING_PHYSICAL)
using X_adv = typename AdvectionPhysicalDomain<CzarnyToCartesian<X, Y, R, Theta>>::X_adv;
using Y_adv = typename AdvectionPhysicalDomain<CzarnyToCartesian<X, Y, R, Theta>>::Y_adv;

#elif defined(CZARNY_MAPPING_PSEUDO_CARTESIAN)
using X_adv = typename AdvectionPseudoCartesianDomain<CzarnyToCartesian<X, Y, R, Theta>>::X_adv;
using Y_adv = typename AdvectionPseudoCartesianDomain<CzarnyToCartesian<X, Y, R, Theta>>::Y_adv;

#elif defined(DISCRETE_MAPPING_PSEUDO_CARTESIAN)
using X_adv = typename AdvectionPseudoCartesianDomain<
        DiscreteToCartesian<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>>::X_adv;
using Y_adv = typename AdvectionPseudoCartesianDomain<
        DiscreteToCartesian<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>>::Y_adv;
#endif

} //end namespace

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



    // Parameters of the grid. ---------------------------------------------------------------
    double const rmin = PCpp_double(conf_voicexx, ".Mesh.r_min");
    double const rmax = PCpp_double(conf_voicexx, ".Mesh.r_max");
    int const Nr = PCpp_int(conf_voicexx, ".Mesh.r_size");
    int const Nt = PCpp_int(conf_voicexx, ".Mesh.p_size");
    double const dt = PCpp_double(conf_voicexx, ".Mesh.time_step");
    double const final_time = PCpp_double(conf_voicexx, ".Mesh.final_time");
    bool const save_curves = PCpp_bool(conf_voicexx, ".Mesh.save_curves");
    bool const save_feet = PCpp_bool(conf_voicexx, ".Mesh.save_feet");
    PC_tree_destroy(&conf_voicexx);

    if (save_curves or save_feet) {
        fs::create_directory("output");
    }
    if (save_curves) {
        fs::create_directory("output/curves");
    }

    std::cout << "TESTS ON THE ADVECTION OPERATOR "
              << "FOR [rmin, rmax] = [" << rmin << ", " << rmax << "], "
              << "WITH NrxNt = " << Nr << "x" << Nt << " AND dt = " << dt << ": " << std::endl;

    // BUILD GRIDS ------------------------------------------------------------------------------
    // Grid creation of space. ------------------------------------------------------------------
    CoordR const r_min(rmin);
    CoordR const r_max(rmax);
    IdxStepR const r_size(Nr);

    CoordTheta const p_min(0.0);
    CoordTheta const p_max(2.0 * M_PI);
    IdxStepTheta const p_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordTheta> p_knots(p_size + 1);

    std::ofstream file("r_knots.txt");
    r_knots[0] = r_min;
    file << 0 << "  " << double(r_knots[0]) << std::endl;
    for (int i(1); i < r_size; ++i) {
        r_knots[i] = CoordR(rmin + i * dr);

        file << i << "  " << double(r_knots[i]) << std::endl;
    }
    r_knots[r_size] = r_max;
    file << int(r_size) << "  " << double(r_knots[r_size]) << std::endl;
    file.close();

    for (int i(0); i < p_size + 1; ++i) {
        p_knots[i] = CoordTheta(p_min + i * dp);
    }

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(p_knots);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR const interpolation_idx_range_R(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta const interpolation_idx_range_P(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const grid(interpolation_idx_range_R, interpolation_idx_range_P);



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
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(r_min);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(r_max);

    SplineRThetaEvaluatorConstBound spline_evaluator_extrapol(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());


    std::string key;

    // SELECTION OF THE MAPPING AND THE ADVECTION DOMAIN.
#if defined(CIRCULAR_MAPPING_PHYSICAL)
    CircularToCartesian<X, Y, R, Theta> analytical_mapping;
    CircularToCartesian<X, Y, R, Theta> mapping;
    AdvectionPhysicalDomain advection_idx_range(analytical_mapping);
    std::string const mapping_name = "CIRCULAR";
    std::string const idx_range_name = "PHYSICAL";
    key += "circular_physical";
#else

    double const czarny_e = 0.3;
    double const czarny_epsilon = 1.4;

#if defined(CZARNY_MAPPING_PHYSICAL)
    CzarnyToCartesian<X, Y, R, Theta> analytical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<X, Y, R, Theta> mapping(czarny_e, czarny_epsilon);
    AdvectionPhysicalDomain advection_idx_range(analytical_mapping);
    std::string const mapping_name = "CZARNY";
    std::string const idx_range_name = "PHYSICAL";
    key += "czarny_physical";

#elif defined(CZARNY_MAPPING_PSEUDO_CARTESIAN)
    CzarnyToCartesian<X, Y, R, Theta> analytical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<X, Y, R, Theta> mapping(czarny_e, czarny_epsilon);
    AdvectionPseudoCartesianDomain advection_idx_range(mapping);
    std::string const mapping_name = "CZARNY";
    std::string const idx_range_name = "PSEUDO CARTESIAN";
    key += "czarny_pseudo_cartesian";

#elif defined(DISCRETE_MAPPING_PSEUDO_CARTESIAN)
    CzarnyToCartesian<X, Y, R, Theta> analytical_mapping(czarny_e, czarny_epsilon);
    DiscreteToCartesian mapping
            = DiscreteToCartesian<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>::
                    analytical_to_discrete(analytical_mapping, builder, spline_evaluator_extrapol);
    AdvectionPseudoCartesianDomain advection_idx_range(mapping);
    std::string const mapping_name = "DISCRETE";
    std::string const idx_range_name = "PSEUDO CARTESIAN";
    key += "discrete_pseudo_cartesian";
#endif
#endif

    key += "-";

    // SELECTION OF THE TIME INTEGRATION METHOD.
#if defined(EULER_METHOD)
    Euler<FieldMemRTheta<CoordRTheta>, DVectorFieldMemRTheta<X_adv, Y_adv>> time_stepper(grid);
    std::string const method_name = "EULER";
    key += "euler";

#elif defined(CRANK_NICOLSON_METHOD)
    double const epsilon_CN = 1e-8;
    CrankNicolson<FieldMemRTheta<CoordRTheta>, DVectorFieldMemRTheta<X_adv, Y_adv>>
            time_stepper(grid, 20, epsilon_CN);
    std::string const method_name = "CRANK NICOLSON";
    key += "crank_nicolson";

#elif defined(RK3_METHOD)
    RK3<FieldMemRTheta<CoordRTheta>, DVectorFieldMemRTheta<X_adv, Y_adv>> time_stepper(grid);
    std::string const method_name = "RK3";
    key += "rk3";

#elif defined(RK4_METHOD)
    RK4<FieldMemRTheta<CoordRTheta>, DVectorFieldMemRTheta<X_adv, Y_adv>> time_stepper(grid);
    std::string const method_name = "RK4";
    key += "rk4";
#endif

    key += "-";

    // SELECTION OF THE SIMULATION.
#if defined(TRANSLATION_SIMULATION)
    TranslationSimulation simulation(mapping, rmin, rmax);
    std::string const simu_type = "TRANSLATION";
    key += "Translation";

#elif defined(ROTATION_SIMULATION)
    RotationSimulation simulation(mapping, rmin, rmax);
    std::string const simu_type = "ROTATION";
    key += "Rotation";

#elif defined(DECENTRED_ROTATION_SIMULATION)
    DecentredRotationSimulation simulation(mapping);
    std::string const simu_type = "DECENTRED ROTATION";
    key += "Decentred_rotation";
#endif

    std::string output_folder = key + "_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }

    std::cout << mapping_name << " MAPPING - " << idx_range_name << " DOMAIN - " << method_name
              << " - " << simu_type << " : " << std::endl;
    simulate(
            mapping,
            analytical_mapping,
            grid,
            time_stepper,
            advection_idx_range,
            simulation,
            interpolator,
            builder,
            spline_evaluator_extrapol,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);
}
