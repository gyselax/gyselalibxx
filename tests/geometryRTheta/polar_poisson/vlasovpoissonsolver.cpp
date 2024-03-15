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

#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polarpoissonsolver.hpp"
#include "test_cases.hpp"
#include "vlasovpoissonsolver.hpp"



using PoissonSolver = PolarSplineFEMPoissonSolver;

#if defined(CIRCULAR_MAPPING)
using Mapping = CircularToCartesian<RDimX, RDimY, RDimR, RDimP>;
#elif defined(CZARNY_MAPPING)
using Mapping = CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP>;
#endif
using DiscreteMapping = DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>;



using LHSFunction = CartesianSolution<Mapping>;
using RHSFunction = ManufacturedPoissonTest<LHSFunction>;


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
    IVectR const r_size(PCpp_int(conf_voicexx, ".Mesh.r_size"));

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(PCpp_int(conf_voicexx, ".Mesh.p_size"));

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordP> p_knots(p_size + 1);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);
    for (int i(0); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    for (int i(0); i < p_size + 1; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<IDimR>(SplineInterpPointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(SplineInterpPointsP::get_sampling());

    IDomainR const interpolation_domain_R(SplineInterpPointsR::get_domain());
    IDomainP const interpolation_domain_P(SplineInterpPointsP::get_domain());
    IDomainRP const grid(interpolation_domain_R, interpolation_domain_P);

    SplineRBuilder const r_builder(interpolation_domain_R);
    SplinePBuilder const p_builder(interpolation_domain_P);
    SplineRPBuilder const builder(grid);

#if defined(CIRCULAR_MAPPING)
    const Mapping mapping;
#elif defined(CZARNY_MAPPING)
    const Mapping mapping(0.3, 1.4);
#endif
    SplineEvaluator2D<BSplinesR, BSplinesP> evaluator(
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>);
    DiscreteMapping const discrete_mapping
            = DiscreteMapping::analytical_to_discrete(mapping, builder, evaluator);

    ddc::init_discrete_space<PolarBSplinesRP>(discrete_mapping, r_builder, p_builder);

    auto const dom_bsplinesRP = builder.spline_domain();

    DFieldRP coeff_alpha(grid);
    DFieldRP coeff_beta(grid);
    DFieldRP x(grid);
    DFieldRP y(grid);

    ddc::for_each(grid, [&](IndexRP const irp) {
        coeff_alpha(irp)
                = std::exp(-std::tanh((ddc::coordinate(ddc::select<IDimR>(irp)) - 0.7) / 0.05));
        coeff_beta(irp) = 1.0 / coeff_alpha(irp);
        const ddc::Coordinate<RDimR, RDimP> coord(ddc::coordinate(irp));
        auto const cartesian_coord = mapping(coord);
        x(irp) = ddc::get<RDimX>(cartesian_coord);
        y(irp) = ddc::get<RDimY>(cartesian_coord);
    });

    Spline2D coeff_alpha_spline(dom_bsplinesRP);
    Spline2D coeff_beta_spline(dom_bsplinesRP);

    builder(coeff_alpha_spline, coeff_alpha);
    builder(coeff_beta_spline, coeff_beta);

    Spline2D x_spline_representation(dom_bsplinesRP);
    Spline2D y_spline_representation(dom_bsplinesRP);

    builder(x_spline_representation, x);
    builder(y_spline_representation, y);

    end_time = std::chrono::system_clock::now();
    std::cout << "Setup time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;
    start_time = std::chrono::system_clock::now();

    PoissonSolver const poisson_solver(coeff_alpha_spline, coeff_beta_spline, discrete_mapping);
    VlasovPoissonSolver<Mapping> const
            vlasov_poisson_solver(mapping, builder, evaluator, poisson_solver);


    end_time = std::chrono::system_clock::now();
    std::cout << "Poisson initialisation time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;


    LHSFunction lhs(mapping);
    RHSFunction rhs(mapping);
    FieldRP<CoordRP> coords(grid);
    DFieldRP result(grid);
    VectorField<double, IDomainRP, NDTag<RDimX, RDimY>> electric_field(grid);

    ddc::for_each(grid, [&](IndexRP const irp) { coords(irp) = CoordRP(ddc::coordinate(irp)); });


    DFieldRP rhs_vals(grid);
    ddc::for_each(grid, [&](IndexRP const irp) { rhs_vals(irp) = rhs(coords(irp)); });

    vlasov_poisson_solver(result.span_view(), electric_field, rhs_vals);


    std::cout << "Solver time : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                         .count()
              << "ms" << std::endl;


    double max_err_x = 0.0;
    double max_err_y = 0.0;
    ddc::for_each(grid, [&](IndexRP const irp) {
        const double diff_x
                = fabs(ddcHelper::get<RDimX>(electric_field)(irp) + lhs.derivative_x(coords(irp)));
        const double diff_y
                = fabs(ddcHelper::get<RDimY>(electric_field)(irp) + lhs.derivative_y(coords(irp)));
        max_err_x = max_err_x > diff_x ? max_err_x : diff_x;
        max_err_y = max_err_y > diff_y ? max_err_y : diff_y;
    });
    std::cout << "Max error Ex : " << max_err_x << std::endl;
    std::cout << "Max error Ey : " << max_err_y << std::endl;


    double max_err = 0.0;
    ddc::for_each(grid, [&](IndexRP const irp) {
        const double err = fabs(result(irp) - lhs(coords(irp)));
        max_err = max_err > err ? max_err : err;
    });
    std::cout << "Max error function : " << max_err << std::endl << std::endl;

    PC_tree_destroy(&conf_voicexx);
    return 0;
}
