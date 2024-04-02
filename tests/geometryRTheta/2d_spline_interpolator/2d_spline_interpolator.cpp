/*
 * Testing the "src/interpolation/spline_interpolator_2d_rp.hpp" file.
 */

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "spline_interpolator_2d_rp.hpp"


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
 * @see SplineInterpolatorRP
 */
template <class Function>
void Interpolation_on_random_coord(
        IDomainRP& grid,
        bool& On_the_nodes,
        Function& exact_function,
        double const TOL)
{
    int const Nr = ddc::discrete_space<BSplinesR>().ncells();
    int const Nt = ddc::discrete_space<BSplinesP>().ncells();


    // Evaluation of the function on the grid. -----------------------------------------------
    DFieldRP function_evaluated(grid);
    ddc::for_each(grid, [&](IndexRP const irp) {
        CoordRP coord(coordinate(ddc::select<IDimR>(irp)), coordinate(ddc::select<IDimP>(irp)));
        function_evaluated(irp) = exact_function(coord);
    });


    // Build the decomposition of the function on Bsplines. ----------------------------------
    SplineRPBuilder const builder(grid);


    // Build a "random" grid to test the interpolator. ---------------------------------------
    IDomainRP random_grid(grid);
    FieldRP<CoordRP> random_coords(random_grid);

    int number_r(0);
    int number_p(0);
    double random_factor_r;
    double random_factor_p;
    ddc::for_each(random_grid, [&](IndexRP const irp) {
        CoordR coord_r(coordinate(ddc::select<IDimR>(irp)));
        CoordP coord_p(coordinate(ddc::select<IDimP>(irp)));

        if (On_the_nodes) {
            random_coords(irp) = CoordRP(coord_r, coord_p);
        } else {
            unsigned index_r = ddc::select<IDimR>(irp).uid();
            unsigned index_p = ddc::select<IDimP>(irp).uid();
            number_r = (301 * index_r + 3) % (Nr * Nt);
            number_p = (103 * index_p + 2) % (Nr * Nt);
            // Pseudo-random number between 0 and 1 generated :
            random_factor_r = double(number_r) / Nr / Nt * (1 - double(number_p) / Nr / Nt);
            random_factor_p = double(number_p) / Nr / Nt * (1 - double(number_r) / Nr / Nt);

            IndexR ir(ddc::select<IDimR>(irp));
            IndexP ip(ddc::select<IDimP>(irp));

            auto r_domain = ddc::get_domain<IDimR>(random_coords);
            IndexR ir_max(r_domain.back());
            CoordR delta_coord_r;
            if (ir + 1 <= ir_max) {
                delta_coord_r = CoordR(ddc::coordinate(ir + 1) - ddc::coordinate(ir));
            } else {
                delta_coord_r = CoordR(0.);
            }

            auto theta_domain = ddc::get_domain<IDimP>(random_coords);
            IndexP ip_min(theta_domain.front());
            IndexP ip_max(theta_domain.back());
            CoordP delta_coord_p;
            if (ip + 1 <= ip_max) {
                delta_coord_p = CoordP(ddc::coordinate(ip + 1) - ddc::coordinate(ip));
            } else {
                delta_coord_p = CoordP(ddc::coordinate(ip_min) - ddc::coordinate(ip_min + 1));
            }

            random_coords(irp) = CoordRP(
                    coord_r + delta_coord_r * random_factor_r,
                    CoordP(fmod(coord_p + delta_coord_p * random_factor_p, 2 * M_PI)));
        }
    });


    // Interpolate the function on Bsplines on the "random" grid. ----------------------------
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<RDimP> p_extrapolation_rule;
    SplineRPEvaluatorNullBound spline_evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    DSpanRP function_interpolated;
    SplineInterpolatorRP interpolator(builder, spline_evaluator);
    function_interpolated = interpolator(function_evaluated, random_coords);

    // Compare the obtained values with the exact function. ----------------------------------
    double max_err(0.0);
    ddc::for_each(random_coords.domain(), [&](IndexRP const irp) {
        double err = fabs(function_interpolated(irp) - exact_function(random_coords(irp)));
        max_err = max_err > err ? max_err : err;
    });
    std::cout << "   Max absolute error : " << max_err;
    if (max_err < TOL) {
        std::cout << "  PASSED" << std::endl;
    } else {
        std::cout << "  FAILED" << std::endl;
    }
}

/**
 * @brief A class of function on polar coordinates.
 */
class ExactFunction
{
public:
    virtual ~ExactFunction() = default;

    /**
     * @brief Evaluate the function.
     * @param[in] coord Point in polar coordinates.
     * @return The value (double) of the function at the point.
     */
    virtual double operator()(CoordRP const& coord) const = 0;
};


/**
 * @brief A class of r polynomial functions on polar coordinates.
 */
class ExactFunction_r_degree : public ExactFunction
{
private:
    const int m_d;

public:
    ExactFunction_r_degree(int d) : m_d(d) {};
    ~ExactFunction_r_degree() {};

    /**
     * @brief Evaluate the function.
     * @param[in] coord Point in polar coordinates.
     * @return The value (double) of the function at the point.
     */
    double operator()(CoordRP const& coord) const override
    {
        const double r = ddc::get<RDimR>(coord);
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
class ExactFunction_r_cos_theta : public ExactFunction
{
public:
    ExactFunction_r_cos_theta() {};
    ~ExactFunction_r_cos_theta() {};

    /**
     * @brief Evaluate the function.
     * @param[in] coord Point in polar coordinates.
     * @return The value (double) of the function at the point.
     */
    double operator()(CoordRP const& coord) const override
    {
        const double r = ddc::get<RDimR>(coord);
        const double t = ddc::get<RDimP>(coord);
        return r * std::cos(t);
    }
};

/**
 * @brief A class of Gaussian function on polar coordinates.
 */
class ExactFunction_gaussian : public ExactFunction
{
public:
    ExactFunction_gaussian() {};
    ~ExactFunction_gaussian() {};

    /**
     * @brief Evaluate the function.
     * @param[in] coord Point in polar coordinates.
     * @return The value (double) of the function at the point.
     */
    double operator()(CoordRP const& coord) const override
    {
        const double r = ddc::get<RDimR>(coord);
        const double t = ddc::get<RDimP>(coord);

        double x0 = 0.;
        double y0 = 0.;
        double sig_x = 0.25;
        double sig_y = 0.25;
        return std::exp(
                -pow(r * std::cos(t) - x0, 2) / (2 * sig_x * sig_x)
                - pow(r * std::sin(t) - y0, 2) / (2 * sig_y * sig_y));
        ;
    }
};


} // end namespace



namespace fs = std::filesystem;

/**
 * @brief Test the interpolator on bsplines in polar coordinates. Test on the mesh point
 * and not on the mesh points for several test functions.
 *
 * Test the "src/interpolation/spline_interpolator_2d_rp.hpp" file.
 * The degree of bsplines is defined in the geometry.hpp file. Normally, the degree is set
 * at 3 for the bsplines in the r dimension and the theta dimension. The interpolation on
 * mesh points must be exact. The interpolation of r-polynomial functions of degree d <= 3
 * must be exact. Otherwise, the interpolation must have a convergence order of 4.
 *
 * @param[in] The path to the grid_size.yaml file.
 *
 * @see SplineInterpolatorRP
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
    int Nt = PCpp_int(conf_voicexx, ".Mesh.p_size");


    int spline_degree = BSplinesR::degree();

    // Grid creation (uniform grid). ----------------------------------------------------------
    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
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
    for (int i(0); i < p_size + 1; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }


    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<IDimR>(SplineInterpPointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(SplineInterpPointsP::get_sampling());

    IDomainR interpolation_domain_R(SplineInterpPointsR::get_domain());
    IDomainP interpolation_domain_P(SplineInterpPointsP::get_domain());
    IDomainRP grid(interpolation_domain_R, interpolation_domain_P);


    // TESTS ON THE DISCRETE SPACE ===========================================================
    std::cout << "Test interpolator 2D with Bsplines of degree " << spline_degree
              << " on r and degree " << BSplinesP::degree() << " on theta"
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
