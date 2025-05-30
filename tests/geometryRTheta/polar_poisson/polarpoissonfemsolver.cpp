// SPDX-License-Identifier: MIT
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "geometry.hpp"
#include "mesh_builder.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polarpoissonlikesolver.hpp"
#include "test_cases.hpp"


using PoissonSolver = PolarSplineFEMPoissonLikeSolver<
        GridR,
        GridTheta,
        PolarBSplinesRTheta,
        SplineRThetaEvaluatorNullBound>;

#if defined(CIRCULAR_MAPPING)
using Mapping = CircularToCartesian<R, Theta, X, Y>;
#elif defined(CZARNY_MAPPING)
using Mapping = CzarnyToCartesian<R, Theta, X, Y>;
#endif
using DiscreteMappingBuilder
        = DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorNullBound>;

using DiscreteMappingBuilder_host = DiscreteToCartesianBuilder<
        X,
        Y,
        SplineRThetaBuilder_host,
        SplineRThetaEvaluatorNullBound_host>;

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

    std::chrono::time_point<std::chrono::system_clock> start_time
            = std::chrono::system_clock::now();
    std::chrono::time_point<std::chrono::system_clock> end_time;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_ncells(PCpp_int(conf_gyselalibxx, ".SplineMesh.r_ncells"));

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(PCpp_int(conf_gyselalibxx, ".SplineMesh.theta_ncells"));

    std::vector<CoordR> r_knots = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> theta_knots
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesR>(r_knots);

    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR interpolation_idx_range_r(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_theta(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_r, interpolation_idx_range_theta);

    SplineRThetaBuilder const builder(grid);
    SplineRThetaBuilder_host const builder_host(grid);

    double major_radius = 6.1;
    double vertical_offset = 0.3;
#if defined(CIRCULAR_MAPPING)
    const Mapping mapping(major_radius, vertical_offset);
#elif defined(CZARNY_MAPPING)
    const Mapping mapping(0.3, 1.4, major_radius, vertical_offset);
#endif

    ddc::NullExtrapolationRule bv_r_min;
    ddc::NullExtrapolationRule bv_r_max;
    ddc::PeriodicExtrapolationRule<Theta> bv_theta_min;
    ddc::PeriodicExtrapolationRule<Theta> bv_theta_max;
    SplineRThetaEvaluatorNullBound evaluator(bv_r_min, bv_r_max, bv_theta_min, bv_theta_max);
    SplineRThetaEvaluatorNullBound_host
            evaluator_host(bv_r_min, bv_r_max, bv_theta_min, bv_theta_max);

    DiscreteMappingBuilder_host const discrete_mapping_builder_host(
            Kokkos::DefaultHostExecutionSpace(),
            mapping,
            builder_host,
            evaluator_host);


    DiscreteMappingBuilder const
            discrete_mapping_builder(Kokkos::DefaultExecutionSpace(), mapping, builder, evaluator);
    DiscreteToCartesian const discrete_mapping = discrete_mapping_builder();
    DiscreteToCartesian const discrete_mapping_host = discrete_mapping_builder_host();

    ddc::init_discrete_space<PolarBSplinesRTheta>(discrete_mapping_host);

    IdxRangeBSRTheta idx_range_bsplinesRTheta = get_spline_idx_range(builder);

    DFieldMemRTheta coeff_alpha_alloc(grid); // values of the coefficient alpha
    DFieldMemRTheta coeff_beta_alloc(grid);
    DFieldMemRTheta x_alloc(grid);
    DFieldMemRTheta y_alloc(grid);

    DFieldRTheta coeff_alpha = get_field(coeff_alpha_alloc); // values of the coefficient alpha
    DFieldRTheta coeff_beta = get_field(coeff_beta_alloc);
    DFieldRTheta x = get_field(x_alloc);
    DFieldRTheta y = get_field(y_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            grid,
            KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                coeff_alpha(irtheta) = Kokkos::exp(
                        -Kokkos::tanh((ddc::coordinate(ddc::select<GridR>(irtheta)) - 0.7) / 0.05));
                coeff_beta(irtheta) = 1.0 / coeff_alpha(irtheta);
                Coord<R, Theta>
                        coord(ddc::coordinate(ddc::select<GridR>(irtheta)),
                              ddc::coordinate(ddc::select<GridTheta>(irtheta)));
                Coord<X, Y> cartesian_coord = mapping(coord);
                x(irtheta) = ddc::get<X>(cartesian_coord);
                y(irtheta) = ddc::get<Y>(cartesian_coord);
            });

    Spline2DMem coeff_alpha_spline(idx_range_bsplinesRTheta);
    Spline2DMem coeff_beta_spline(idx_range_bsplinesRTheta);

    builder(get_field(coeff_alpha_spline),
            get_const_field(coeff_alpha)); // coeff_alpha_spline are the coefficients
    // of the spline representation of the values given by coeff_alpha.
    builder(get_field(coeff_beta_spline), get_const_field(coeff_beta));

    Spline2DMem x_spline_representation(idx_range_bsplinesRTheta);
    Spline2DMem y_spline_representation(idx_range_bsplinesRTheta);

    builder(get_field(x_spline_representation), get_const_field(x));
    builder(get_field(y_spline_representation), get_const_field(y));

    end_time = std::chrono::system_clock::now();
    std::cout << "Setup time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;
    start_time = std::chrono::system_clock::now();

    PoissonSolver
            solver(get_const_field(coeff_alpha_spline),
                   get_const_field(coeff_beta_spline),
                   discrete_mapping,
                   evaluator);

    end_time = std::chrono::system_clock::now();
    std::cout << "Poisson initialisation time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;

    LHSFunction lhs(mapping, major_radius, vertical_offset);
    RHSFunction rhs(mapping);
    host_t<FieldMemRTheta<CoordRTheta>> coords(grid);
    DFieldMemRTheta result_alloc(grid);
    DField<IdxRangeRTheta> result = get_field(result_alloc);

    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        coords(irtheta) = CoordRTheta(
                ddc::coordinate(ddc::select<GridR>(irtheta)),
                ddc::coordinate(ddc::select<GridTheta>(irtheta)));
    });
    if (discrete_rhs) {
        // Build the spline approximation of the rhs
        Spline2DMem rhs_spline(idx_range_bsplinesRTheta);
        host_t<DFieldMemRTheta> rhs_vals_host(grid);

        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            rhs_vals_host(irtheta) = rhs(coords(irtheta));
        });
        auto rhs_vals = ddc::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(rhs_vals_host));

        builder(get_field(rhs_spline), get_const_field(rhs_vals));
        ConstSpline2D rhs_spline_field = get_const_field(rhs_spline);
        start_time = std::chrono::system_clock::now();
        solver(
                KOKKOS_LAMBDA(CoordRTheta const& coord) {
                    return evaluator(coord, rhs_spline_field);
                },
                get_field(result));
        end_time = std::chrono::system_clock::now();
    } else {
        start_time = std::chrono::system_clock::now();
        solver(rhs, get_field(result));
        end_time = std::chrono::system_clock::now();
    }
    auto result_alloc_host = ddc::create_mirror_view_and_copy(result);
    host_t<DFieldRTheta> result_host = get_field(result_alloc_host);
    std::cout << "Solver time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;

    double max_err = 0.0;
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        const double err = result_host(irtheta) - lhs(coords(irtheta));
        if (err > 0) {
            max_err = max_err > err ? max_err : err;
        } else {
            max_err = max_err > -err ? max_err : -err;
        }
    });
    std::cout << "Max error : " << max_err << std::endl;

    PC_tree_destroy(&conf_gyselalibxx);
    return 0;
}
