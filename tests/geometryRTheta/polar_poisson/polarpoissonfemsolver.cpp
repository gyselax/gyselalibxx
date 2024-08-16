// SPDX-License-Identifier: MIT
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "mesh_builder.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polarpoissonlikesolver.hpp"
#include "test_cases.hpp"

using PoissonSolver = PolarSplineFEMPoissonLikeSolver;

#if defined(CIRCULAR_MAPPING)
using Mapping = CircularToCartesian<X, Y, R, Theta>;
#elif defined(CZARNY_MAPPING)
using Mapping = CzarnyToCartesian<X, Y, R, Theta>;
#endif
using DiscreteMapping
        = DiscreteToCartesian<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorNullBound>;

#if defined(CURVILINEAR_SOLUTION)
using LHSFunction = CurvilinearSolution<Mapping>;
#elif defined(CARTESIAN_SOLUTION)
using LHSFunction = CartesianSolution<Mapping>;
#endif
using RHSFunction = ManufacturedPoissonTest<LHSFunction>;

constexpr bool discrete_rhs = false;

namespace fs = std::filesystem;

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

    std::chrono::time_point<std::chrono::system_clock> start_time
            = std::chrono::system_clock::now();
    std::chrono::time_point<std::chrono::system_clock> end_time;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_ncells(PCpp_int(conf_voicexx, ".SplineMesh.r_ncells"));

    CoordTheta const p_min(0.0);
    CoordTheta const p_max(2.0 * M_PI);
    IdxStepTheta const p_ncells(PCpp_int(conf_voicexx, ".SplineMesh.p_ncells"));

    std::vector<CoordR> r_knots = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> p_knots = build_uniform_break_points(p_min, p_max, p_ncells);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesR>(r_knots);

    ddc::init_discrete_space<BSplinesTheta>(p_knots);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR interpolation_idx_range_R(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_P(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_R, interpolation_idx_range_P);

    SplineRThetaBuilder const builder(grid);

#if defined(CIRCULAR_MAPPING)
    const Mapping mapping;
#elif defined(CZARNY_MAPPING)
    const Mapping mapping(0.3, 1.4);
#endif
    ddc::NullExtrapolationRule bv_r_min;
    ddc::NullExtrapolationRule bv_r_max;
    ddc::PeriodicExtrapolationRule<Theta> bv_p_min;
    ddc::PeriodicExtrapolationRule<Theta> bv_p_max;
    SplineRThetaEvaluatorNullBound evaluator(bv_r_min, bv_r_max, bv_p_min, bv_p_max);
    DiscreteMapping const discrete_mapping
            = DiscreteMapping::analytical_to_discrete(mapping, builder, evaluator);

    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping);

    auto dom_bsplinesRTheta = get_spline_idx_range(builder);

    DFieldMemRTheta coeff_alpha(grid); // values of the coefficent alpha
    DFieldMemRTheta coeff_beta(grid);
    DFieldMemRTheta x(grid);
    DFieldMemRTheta y(grid);

    ddc::for_each(grid, [&](IdxRTheta const irp) {
        coeff_alpha(irp)
                = std::exp(-std::tanh((ddc::coordinate(ddc::select<GridR>(irp)) - 0.7) / 0.05));
        coeff_beta(irp) = 1.0 / coeff_alpha(irp);
        Coord<R, Theta>
                coord(ddc::coordinate(ddc::select<GridR>(irp)),
                      ddc::coordinate(ddc::select<GridTheta>(irp)));
        auto cartesian_coord = mapping(coord);
        x(irp) = ddc::get<X>(cartesian_coord);
        y(irp) = ddc::get<Y>(cartesian_coord);
    });

    Spline2D coeff_alpha_spline(dom_bsplinesRTheta);
    Spline2D coeff_beta_spline(dom_bsplinesRTheta);

    builder(get_field(coeff_alpha_spline),
            get_const_field(coeff_alpha)); // coeff_alpha_spline are the coefficients
    // of the spline representation of the values given by coeff_alpha.
    builder(get_field(coeff_beta_spline), get_const_field(coeff_beta));

    Spline2D x_spline_representation(dom_bsplinesRTheta);
    Spline2D y_spline_representation(dom_bsplinesRTheta);

    builder(get_field(x_spline_representation), get_const_field(x));
    builder(get_field(y_spline_representation), get_const_field(y));

    end_time = std::chrono::system_clock::now();
    std::cout << "Setup time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;
    start_time = std::chrono::system_clock::now();

    PoissonSolver solver(coeff_alpha_spline, coeff_beta_spline, discrete_mapping);

    end_time = std::chrono::system_clock::now();
    std::cout << "Poisson initialisation time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;

    LHSFunction lhs(mapping);
    RHSFunction rhs(mapping);
    FieldMemRTheta<CoordRTheta> coords(grid);
    DFieldMemRTheta result(grid);
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        coords(irp) = CoordRTheta(
                ddc::coordinate(ddc::select<GridR>(irp)),
                ddc::coordinate(ddc::select<GridTheta>(irp)));
    });
    if (discrete_rhs) {
        // Build the spline approximation of the rhs

        Spline2D rhs_spline(dom_bsplinesRTheta);
        DFieldMemRTheta rhs_vals(grid);
        ddc::for_each(grid, [&](IdxRTheta const irp) { rhs_vals(irp) = rhs(coords(irp)); });
        builder(get_field(rhs_spline), get_const_field(rhs_vals));



        start_time = std::chrono::system_clock::now();
        ddc::NullExtrapolationRule r_extrapolation_rule;
        ddc::PeriodicExtrapolationRule<Theta> p_extrapolation_rule;
        SplineRThetaEvaluatorNullBound
                eval(r_extrapolation_rule,
                     r_extrapolation_rule,
                     p_extrapolation_rule,
                     p_extrapolation_rule);
        solver([&](CoordRTheta const& coord) { return eval(coord, get_const_field(rhs_spline)); },
               get_const_field(coords),
               get_field(result));
        end_time = std::chrono::system_clock::now();
    } else {
        start_time = std::chrono::system_clock::now();
        solver(rhs, get_const_field(coords), get_field(result));
        end_time = std::chrono::system_clock::now();
    }
    std::cout << "Solver time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;

    double max_err = 0.0;
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        const double err = result(irp) - lhs(coords(irp));
        if (err > 0) {
            max_err = max_err > err ? max_err : err;
        } else {
            max_err = max_err > -err ? max_err : -err;
        }
    });
    std::cout << "Max error : " << max_err << std::endl;

    PC_tree_destroy(&conf_voicexx);
    return 0;
}
