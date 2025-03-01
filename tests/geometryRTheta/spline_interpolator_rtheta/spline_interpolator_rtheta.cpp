/*
 * Testing the "src/interpolation/spline_interpolator_rtheta.hpp" file.
 */

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "spline_interpolator_rtheta.hpp"


/*
 * Remark: instead of working with pseudo-random interpolation grid,
 *         using #include <random> could be implemented.
 *         see: https://en.cppreference.com/w/cpp/numeric/random
 */


namespace {


/**
 * @brief Interpolate a function on bsplines in polar coordinates and
 * print the maximum absolute error of the interpolation.
 *
 * @param[in] grid The polar discretisation grid on which we are working.
 * @param[in] On_the_nodes
 * 				If true, interpolate on the mesh points.
 * 				If False, interpolate on another grid with points located
 * 				at pseudo-random place in each cell.
 * @param[in] exact_function The exact function that we want to interpolate.
 * @param[in] TOL A tolerance threshold above which the error is considered
 * 				to be sufficiently small.
 *
 * @see SplineInterpolatorRTheta
 */
template <class Function>
void Interpolation_on_random_coord(
        IdxRangeRTheta& grid,
        bool& On_the_nodes,
        Function& exact_function,
        double const TOL)
{
    int const Nr = ddc::discrete_space<BSplinesR>().ncells();
    int const Nt = ddc::discrete_space<BSplinesTheta>().ncells();


    // Evaluation of the function on the grid. -----------------------------------------------
    DFieldMemRTheta function_evaluated_alloc(grid);
    DFieldRTheta function_evaluated = get_field(function_evaluated_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            grid,
            KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                CoordRTheta coord(ddc::coordinate(irtheta));
                function_evaluated(irtheta) = exact_function(coord);
            });


    // Build the decomposition of the function on Bsplines. ----------------------------------
    SplineRThetaBuilder const builder(grid);


    // Build a "random" grid to test the interpolator. ---------------------------------------
    FieldMemRTheta<CoordRTheta> random_coords_alloc(grid);
    FieldRTheta<CoordRTheta> random_coords = get_field(random_coords_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            grid,
            KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                CoordR coord_r(coordinate(ddc::select<GridR>(irtheta)));
                CoordTheta coord_theta(coordinate(ddc::select<GridTheta>(irtheta)));

                if (On_the_nodes) {
                    random_coords(irtheta) = CoordRTheta(coord_r, coord_theta);
                } else {
                    IdxRangeR r_idx_range(grid);
                    IdxRangeTheta theta_idx_range(grid);

                    IdxR ir_min(r_idx_range.front());
                    IdxTheta itheta_min(theta_idx_range.front());

                    IdxR ir_max(r_idx_range.back());
                    IdxTheta itheta_max(theta_idx_range.back());

                    IdxR ir(irtheta);
                    IdxTheta itheta(irtheta);

                    unsigned index_r = ir - ir_min;
                    unsigned index_theta = itheta - itheta_min;
                    int const number_r = (301 * index_r + 3) % (Nr * Nt);
                    int const number_theta = (103 * index_theta + 2) % (Nr * Nt);
                    // Pseudo-random number between 0 and 1 generated :
                    double const random_factor_r
                            = double(number_r) / Nr / Nt * (1 - double(number_theta) / Nr / Nt);
                    double const random_factor_theta
                            = double(number_theta) / Nr / Nt * (1 - double(number_r) / Nr / Nt);

                    CoordR delta_coord_r;
                    if (ir + 1 <= ir_max) {
                        delta_coord_r = CoordR(ddc::coordinate(ir + 1) - ddc::coordinate(ir));
                    } else {
                        delta_coord_r = CoordR(0.);
                    }

                    CoordTheta delta_coord_theta;
                    if (itheta + 1 <= itheta_max) {
                        delta_coord_theta = CoordTheta(ddc::coordinate(itheta + 1) - ddc::coordinate(itheta));
                    } else {
                        delta_coord_theta
                                = CoordTheta(ddc::coordinate(itheta_min) - ddc::coordinate(itheta_min + 1));
                    }

                    random_coords(irtheta) = CoordRTheta(
                            coord_r + delta_coord_r * random_factor_r,
                            CoordTheta(fmod(coord_theta + delta_coord_theta * random_factor_theta, 2 * M_PI)));
                }
            });


    // Interpolate the function on Bsplines on the "random" grid. ----------------------------
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

    SplineInterpolatorRTheta interpolator(builder, spline_evaluator);
    interpolator(get_field(function_evaluated), get_const_field(random_coords));

    // Compare the obtained values with the exact function. ----------------------------------
    double max_err = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            grid,
            0.0,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                return Kokkos::fabs(function_evaluated(irtheta) - exact_function(random_coords(irtheta)));
            });

    std::cout << "   Max absolute error : " << max_err;
    if (max_err < TOL) {
        std::cout << "  PASSED" << std::endl;
    } else {
        std::cout << "  FAILED" << std::endl;
    }
}


/**
 * @brief A class of r polynomial functions on polar coordinates.
 */
class ExactFunction_r_degree
{
private:
    const int m_d;

public:
    explicit KOKKOS_FUNCTION ExactFunction_r_degree(int d) : m_d(d) {}

    /**
     * @brief Evaluate the function.
     * @param[in] coord Point in polar coordinates.
     * @return The value (double) of the function at the point.
     */
    KOKKOS_FUNCTION double operator()(CoordRTheta const& coord) const
    {
        const double r = ddc::get<R>(coord);
        double val = 1.0;
        for (int i(0); i < m_d; i++)
            val *= r;
        return val;
    }
};

/**
 * @brief A class of function on polar coordinates such that
 * f(r, theta) = r *cos(theta).
 */
class ExactFunction_r_cos_theta
{
public:
    /**
     * @brief Evaluate the function.
     * @param[in] coord Point in polar coordinates.
     * @return The value (double) of the function at the point.
     */
    KOKKOS_FUNCTION double operator()(CoordRTheta const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);
        return r * Kokkos::cos(theta);
    }
};

/**
 * @brief A class of Gaussian function on polar coordinates.
 */
class ExactFunction_gaussian
{
public:
    /**
     * @brief Evaluate the function.
     * @param[in] coord Point in polar coordinates.
     * @return The value (double) of the function at the point.
     */
    KOKKOS_FUNCTION double operator()(CoordRTheta const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double theta = ddc::get<Theta>(coord);

        double x0 = 0.;
        double y0 = 0.;
        double sig_x = 0.25;
        double sig_y = 0.25;
        return Kokkos::exp(
                -ipow(r * Kokkos::cos(theta) - x0, 2) / (2 * sig_x * sig_x)
                - ipow(r * Kokkos::sin(theta) - y0, 2) / (2 * sig_y * sig_y));
        ;
    }
};


} // end namespace



namespace fs = std::filesystem;

/**
 * @brief Test the interpolator on bsplines in polar coordinates. Test on the mesh point
 * and not on the mesh points for several test functions.
 *
 * Test the "src/interpolation/spline_interpolator_rtheta.hpp" file.
 * The degree of bsplines is defined in the geometry.hpp file. Normally, the degree is set
 * at 3 for the bsplines in the r dimension and the theta dimension. The interpolation on
 * mesh points must be exact. The interpolation of r-polynomial functions of degree d <= 3
 * must be exact. Otherwise, the interpolation must have a convergence order of 4.
 *
 * @param[in] The path to the grid_size.yaml file.
 *
 * @see SplineInterpolatorRTheta
 * @tag Spline_interpolator_polar_test
 */
int main(int argc, char** argv)
{
    // INITIALISATION OF THE DISCRETE SPACE ==================================================
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
    int Nr = PCpp_int(conf_voicexx, ".Mesh.r_size");
    int Nt = PCpp_int(conf_voicexx, ".Mesh.theta_size");


    int spline_degree = BSplinesR::degree();

    // Grid creation (uniform grid). ----------------------------------------------------------
    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dp((theta_max - theta_min) / theta_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordTheta> theta_knots(theta_size + 1);

    for (int i(0); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    for (int i(0); i < theta_size + 1; ++i) {
        theta_knots[i] = CoordTheta(theta_min + i * dp);
    }


    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR interpolation_idx_range_R(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_Theta(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_R, interpolation_idx_range_Theta);


    // TESTS ON THE DISCRETE SPACE ===========================================================
    std::cout << "Test interpolator 2D with Bsplines of degree " << spline_degree
              << " on r and degree " << BSplinesTheta::degree() << " on theta"
              << " on a grid of " << Nr << " x " << Nt << ":" << std::endl;

    // For interpolation points on the node, the interpolation is supposed to be exact:
    bool On_the_nodes(true);
    std::cout << "- For interpolation points on the node grid: ------------------------"
              << std::endl;

    for (int degree(spline_degree); degree <= spline_degree + 4; degree++) {
        std::cout << "   Test with f(r,th) = r^" << degree << " :      ";
        ExactFunction_r_degree exact_function_Rd(degree);
        Interpolation_on_random_coord<
                ExactFunction_r_degree>(grid, On_the_nodes, exact_function_Rd, 1e-15);
    }

    std::cout << "  -------------------------------------------------------------------"
              << std::endl;

    std::cout << "   Test with f(r,th) = r*cos(th) :";
    ExactFunction_r_cos_theta exact_function_r_cos_theta;
    Interpolation_on_random_coord<
            ExactFunction_r_cos_theta>(grid, On_the_nodes, exact_function_r_cos_theta, 1e-15);

    std::cout << "  -------------------------------------------------------------------"
              << std::endl;


    std::cout << "   Test with f(r,th) = Gaussian  :";
    ExactFunction_gaussian exact_function_gaussian;
    Interpolation_on_random_coord<
            ExactFunction_gaussian>(grid, On_the_nodes, exact_function_gaussian, 1e-15);

    std::cout << "  -------------------------------------------------------------------"
              << std::endl;


    // For interpolation points not on the node, the interpolation is supposed to be exact
    // or order d+1:
    On_the_nodes = false;
    std::cout << "- For interpolation points not on the node grid: --------------------"
              << std::endl;

    std::cout << "   Test with f(r,th) = 1.0 :      ";
    ExactFunction_r_degree exact_function_R0(0);
    Interpolation_on_random_coord<
            ExactFunction_r_degree>(grid, On_the_nodes, exact_function_R0, 5e-15);

    for (int degree(1); degree <= spline_degree; degree++) {
        std::cout << "   Test with f(r,th) = r^" << degree << " :      ";
        ExactFunction_r_degree exact_function_Rd(degree);
        Interpolation_on_random_coord<
                ExactFunction_r_degree>(grid, On_the_nodes, exact_function_Rd, 5e-15);
    }

    for (int degree(spline_degree + 1); degree <= spline_degree + 4; degree++) {
        std::cout << "   Test with f(r,th) = r^" << degree << " :      ";
        ExactFunction_r_degree exact_function_Rd(degree);
        Interpolation_on_random_coord<ExactFunction_r_degree>(
                grid,
                On_the_nodes,
                exact_function_Rd,
                5e-3 * std::pow(16. * 16. / Nr / Nt, spline_degree / 2.));
    }

    std::cout << "  -------------------------------------------------------------------"
              << std::endl;

    std::cout << "   Test with f(r,th) = r*cos(th) :";
    //ExactFunction_r_cos_theta exact_function_r_cos_theta; //(already defined)
    Interpolation_on_random_coord<ExactFunction_r_cos_theta>(
            grid,
            On_the_nodes,
            exact_function_r_cos_theta,
            5e-3 * std::pow(16. * 16. / Nr / Nt, spline_degree / 2.));

    std::cout << "  -------------------------------------------------------------------"
              << std::endl;


    std::cout << "   Test with f(r,th) = Gaussian  :";
    //ExactFunction_gaussian exact_function_gaussian; //(already defined)
    Interpolation_on_random_coord<ExactFunction_gaussian>(
            grid,
            On_the_nodes,
            exact_function_gaussian,
            5e-3 * std::pow(16. * 16. / Nr / Nt, spline_degree / 2.));

    std::cout << "  -------------------------------------------------------------------"
              << std::endl;

    PC_tree_destroy(&conf_voicexx);
}
