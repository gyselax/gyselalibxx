
#include <chrono>
#include <typeinfo>

#include <ddc/ddc.hpp>

#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/math_tools.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include "advection_domain.hpp"
#include "bsl_advection_rp.hpp"
#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "test_cases.hpp"

// ...
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <stdio.h>
// ...

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <crank_nicolson.hpp>
#include <ddc_helper.hpp>
#include <euler.hpp>
#include <rk2.hpp>
#include <rk3.hpp>
#include <rk4.hpp>

#include "advection_simulation_utils.hpp"
#include "itimestepper.hpp"
#include "spline_foot_finder.hpp"



namespace fs = std::filesystem;

namespace {
#if defined(CIRCULAR_MAPPING_PHYSICAL)
using RDimX_adv = typename AdvectionPhysicalDomain<
        CircularToCartesian<RDimX, RDimY, RDimR, RDimP>>::RDimX_adv;
using RDimY_adv = typename AdvectionPhysicalDomain<
        CircularToCartesian<RDimX, RDimY, RDimR, RDimP>>::RDimY_adv;
#elif defined(CZARNY_MAPPING_PHYSICAL)
using RDimX_adv =
        typename AdvectionPhysicalDomain<CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP>>::RDimX_adv;
using RDimY_adv =
        typename AdvectionPhysicalDomain<CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP>>::RDimY_adv;

#elif defined(CZARNY_MAPPING_PSEUDO_CARTESIAN)
using RDimX_adv = typename AdvectionPseudoCartesianDomain<
        CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP>>::RDimX_adv;
using RDimY_adv = typename AdvectionPseudoCartesianDomain<
        CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP>>::RDimY_adv;

#elif defined(DISCRETE_MAPPING_PSEUDO_CARTESIAN)
using RDimX_adv = typename AdvectionPseudoCartesianDomain<
        DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>>::RDimX_adv;
using RDimY_adv = typename AdvectionPseudoCartesianDomain<
        DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>>::RDimY_adv;
#endif

} //end namespace

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
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordP> p_knots(p_size + 1);

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


    std::string key;

    // SELECTION OF THE MAPPING AND THE ADVECTION DOMAIN.
#if defined(CIRCULAR_MAPPING_PHYSICAL)
    CircularToCartesian<RDimX, RDimY, RDimR, RDimP> analytical_mapping;
    CircularToCartesian<RDimX, RDimY, RDimR, RDimP> mapping;
    AdvectionPhysicalDomain advection_domain(analytical_mapping);
    std::string const mapping_name = "CIRCULAR";
    std::string const domain_name = "PHYSICAL";
    key += "circular_physical";
#else

    double const czarny_e = 0.3;
    double const czarny_epsilon = 1.4;

#if defined(CZARNY_MAPPING_PHYSICAL)
    CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP> analytical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP> mapping(czarny_e, czarny_epsilon);
    AdvectionPhysicalDomain advection_domain(analytical_mapping);
    std::string const mapping_name = "CZARNY";
    std::string const domain_name = "PHYSICAL";
    key += "czarny_physical";

#elif defined(CZARNY_MAPPING_PSEUDO_CARTESIAN)
    CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP> analytical_mapping(czarny_e, czarny_epsilon);
    CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP> mapping(czarny_e, czarny_epsilon);
    AdvectionPseudoCartesianDomain advection_domain(mapping);
    std::string const mapping_name = "CZARNY";
    std::string const domain_name = "PSEUDO CARTESIAN";
    key += "czarny_pseudo_cartesian";

#elif defined(DISCRETE_MAPPING_PSEUDO_CARTESIAN)
    CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP> analytical_mapping(czarny_e, czarny_epsilon);
    DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder> mapping
            = DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>::
                    analytical_to_discrete(analytical_mapping, builder, spline_evaluator_extrapol);
    AdvectionPseudoCartesianDomain advection_domain(mapping);
    std::string const mapping_name = "DISCRETE";
    std::string const domain_name = "PSEUDO CARTESIAN";
    key += "discrete_pseudo_cartesian";
#endif
#endif

    key += "-";

    // SELECTION OF THE TIME INTEGRATION METHOD.
#if defined(EULER_METHOD)
    Euler<FieldRP<CoordRP>, VectorDFieldRP<RDimX_adv, RDimY_adv>> time_stepper(grid);
    std::string const method_name = "EULER";
    key += "euler";

#elif defined(CRANK_NICOLSON_METHOD)
    double const epsilon_CN = 1e-8;
    CrankNicolson<FieldRP<CoordRP>, VectorDFieldRP<RDimX_adv, RDimY_adv>>
            time_stepper(grid, 20, epsilon_CN);
    std::string const method_name = "CRANK NICOLSON";
    key += "crank_nicolson";

#elif defined(RK3_METHOD)
    RK3<FieldRP<CoordRP>, VectorDFieldRP<RDimX_adv, RDimY_adv>> time_stepper(grid);
    std::string const method_name = "RK3";
    key += "rk3";

#elif defined(RK4_METHOD)
    RK4<FieldRP<CoordRP>, VectorDFieldRP<RDimX_adv, RDimY_adv>> time_stepper(grid);
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

    std::cout << mapping_name << " MAPPING - " << domain_name << " DOMAIN - " << method_name
              << " - " << simu_type << " : " << std::endl;
    simulate(
            mapping,
            analytical_mapping,
            grid,
            time_stepper,
            advection_domain,
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
