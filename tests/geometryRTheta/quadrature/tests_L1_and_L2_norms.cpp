#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "geometry.hpp"
#include "l_norm_tools.hpp"
#include "mesh_builder.hpp"
#include "quadrature.hpp"
#include "spline_quadrature.hpp"
#include "trapezoid_quadrature.hpp"
#include "volume_quadrature_nd.hpp"



namespace {
/**
 * @brief Check if the L1 norm corresponds to the expected L1 norm.
 *
 * It checks if
 * @f$ |\Vert f \Vert_{\mathcal{L}_1, \text{computed}}
 * - \Vert f \Vert_{\mathcal{L}_1, \text{expected}} | < \varepsilon @f$.
 *
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      The function defined on the polar grid.
 * @param[in] expected_norm
 *      The expected L1 norm of the function.
 * @param[in] TOL
 *      The error tolerance @f$ \varepsilon @f$.
 */
void check_norm_L1(
        host_t<Quadrature<IdxRangeRTheta>> quadrature,
        host_t<DFieldRTheta> function,
        double const expected_norm,
        double const TOL)
{
    double const norm_function = norm_L1(Kokkos::DefaultHostExecutionSpace(), quadrature, function);
    std::cout << "  > error on L1 norm = " << fabs(norm_function - expected_norm) << std::endl;
    EXPECT_NEAR(norm_function, expected_norm, TOL);
}


/**
 * @brief Check if the L2 norm corresponds to the expected L1 norm.
 *
 * It checks if
 * @f$ |\Vert f \Vert_{\mathcal{L}_2, \text{computed}}
 * - \Vert f \Vert_{\mathcal{L}_2, \text{expected}} | < \varepsilon @f$.
 *
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      The function defined on the polar grid.
 * @param[in] expected_norm
 *      The expected L2 norm of the function.
 * @param[in] TOL
 *      The error tolerance @f$ \varepsilon @f$.
 */
void check_norm_L2(
        host_t<Quadrature<IdxRangeRTheta>> quadrature,
        host_t<DFieldRTheta> function,
        double const expected_norm,
        double const TOL)
{
    double const norm_function = norm_L2(Kokkos::DefaultHostExecutionSpace(), quadrature, function);
    std::cout << "  > error on L2 norm = " << fabs(norm_function - expected_norm) << std::endl;
    EXPECT_NEAR(norm_function, expected_norm, TOL);
}


/**
 * @brief Check if the L1 and the L2 norms correspond to the expected norms.
 *
 * Use check_norm_L1 and check_norm_L2.
 *
 * @param[in] quadrature
 *      The quadrature used to compute the integral.
 * @param[in] function
 *      The function defined on the polar grid.
 * @param[in] expected_norms
 *      A 2x1 array with the expected L2 norm of the function.
 * @param[in] TOLs
 *      A 2x1 array with the error tolerance @f$ \varepsilon @f$.
 */
void check_norms(
        host_t<Quadrature<IdxRangeRTheta>> quadrature,
        host_t<DFieldRTheta> function,
        std::array<double, 2> const& expected_norms,
        std::array<double, 2> const& TOLs)
{
    check_norm_L1(quadrature, function, expected_norms[0], TOLs[0]);
    check_norm_L2(quadrature, function, expected_norms[1], TOLs[1]);
}


template <class Mapping>
void launch_tests(
        Mapping const& mapping,
        SplineRThetaBuilder_host const& builder,
        IdxRangeRTheta const& grid,
        std::array<std::array<double, 2>, 5> const& expected_norms,
        std::array<std::array<double, 2>, 5> const& TOLs)
{
    using SplineRBuilder = ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesR,
            GridR,
            ddc::BoundCond::GREVILLE, // boundary at r=0
            ddc::BoundCond::GREVILLE, // boundary at rmax
            ddc::SplineSolver::LAPACK>;
    using SplinePBuilder = ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesTheta,
            GridTheta,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK>;

    SplineRBuilder r_builder(ddc::select<GridR>(grid));
    SplinePBuilder p_builder(ddc::select<GridTheta>(grid));
    // Test spline quadrature: ------------------------------------------------------------------------
    // Instantiate a quadrature with coefficients where we added the Jacobian determinant.
    host_t<DFieldMemRTheta> const quadrature_coeffs = compute_coeffs_on_mapping(
            Kokkos::DefaultHostExecutionSpace(),
            mapping,
            spline_quadrature_coefficients<
                    Kokkos::DefaultHostExecutionSpace>(grid, r_builder, p_builder));
    host_t<Quadrature<IdxRangeRTheta>> quadrature(get_const_field(quadrature_coeffs));

    host_t<DFieldMemRTheta> test(grid);

    // --- TEST 1 -------------------------------------------------------------------------------------
    std::cout << "TEST 1: f(r,theta ) = 1. " << std::endl;
    ddc::for_each(grid, [&](IdxRTheta const irtheta) { test(irtheta) = 1.; });
    check_norms(quadrature, get_field(test), expected_norms[0], TOLs[0]);

    // --- TEST 2 -------------------------------------------------------------------------------------
    std::cout << std::endl << "TEST 2: f(r,theta ) = cos(theta) " << std::endl;
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        double const th = ddc::select<Theta>(ddc::coordinate(irtheta));
        test(irtheta) = std::cos(th);
    });
    check_norms(quadrature, get_field(test), expected_norms[1], TOLs[1]);
    std::cout << "Remark: (r, theta) -> |cos(theta)| not C¹." << std::endl;

    // --- TEST 3 -------------------------------------------------------------------------------------
    std::cout << std::endl << "TEST 3: f(r,theta ) = r " << std::endl;
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        double const r = ddc::select<R>(ddc::coordinate(irtheta));
        test(irtheta) = r;
    });
    check_norms(quadrature, get_field(test), expected_norms[2], TOLs[2]);

    // --- TEST 4 -------------------------------------------------------------------------------------
    std::cout << std::endl << "TEST 4: f(r,theta ) = r cos(theta) " << std::endl;
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        double const r = ddc::select<R>(ddc::coordinate(irtheta));
        double const th = ddc::select<Theta>(ddc::coordinate(irtheta));
        test(irtheta) = r * std::cos(th);
    });
    check_norms(quadrature, get_field(test), expected_norms[3], TOLs[3]);
    std::cout << "Remark: (r, theta) -> |r cos(theta)| not C¹." << std::endl;

    // --- TEST 5 -------------------------------------------------------------------------------------
    std::cout << std::endl << "TEST 5: f(r,theta ) = r⁵ cos(10*theta) " << std::endl;
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        double const r = ddc::select<R>(ddc::coordinate(irtheta));
        double const th = ddc::select<Theta>(ddc::coordinate(irtheta));
        test(irtheta) = std::pow(r, 5) * std::cos(10 * th);
    });
    check_norms(quadrature, get_field(test), expected_norms[4], TOLs[4]);
    std::cout << "Remark: (r, theta) -> |r⁵ cos(10*theta)| not C¹." << std::endl;
}

} // end namespace


/**
 * @brief A class for the Google tests of the SplineQuadratureRTheta.
 */
class SplineQuadrature : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};


TEST_P(SplineQuadrature, TestFunctions)
{
    // Build the grid for the space. ------------------------------------------------------------------
    auto const [Nr, Nt] = GetParam();

    CoordR const r_min(0.);
    CoordR const r_max(1.);
    IdxStepR const r_ncells(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(Nt);

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

    SplineRThetaBuilder_host builder(grid);


    // ------------------------------------------------------------------------------------------------
    std::cout << "CIRCULAR MAPPING ---------------------------------------------------"
              << std::endl;
    const CircularToCartesian<R, Theta, X, Y> mapping_1;
    std::array<std::array<double, 2>, 5> expected_norms;
    expected_norms[0][0] = M_PI;
    expected_norms[0][1] = std::sqrt(M_PI);
    expected_norms[1][0] = 4. / 2.;
    expected_norms[1][1] = std::sqrt(M_PI / 2.);
    expected_norms[2][0] = 2 * M_PI / 3.;
    expected_norms[2][1] = std::sqrt(2 * M_PI / 4.);
    expected_norms[3][0] = 4 / 3.;
    expected_norms[3][1] = std::sqrt(M_PI / 4.);
    expected_norms[4][0] = 4 / 7.;
    expected_norms[4][1] = std::sqrt(M_PI / 12.);

    std::array<std::array<double, 2>, 5> TOLs;
    TOLs[0][0] = 1e-13;
    TOLs[0][1] = 1e-13;
    TOLs[1][0] = 5e-3;
    TOLs[1][1] = 1e-13;
    TOLs[2][0] = 1e-13;
    TOLs[2][1] = 1e-13;
    TOLs[3][0] = 5e-3;
    TOLs[3][1] = 1e-13;
    TOLs[4][0] = 5e-2;
    TOLs[4][1] = 5e-6;

    launch_tests(mapping_1, builder, grid, expected_norms, TOLs);


    std::cout << std::endl
              << "DISCRETE MAPPING ---------------------------------------------------"
              << std::endl;
    ddc::ConstantExtrapolationRule<R, Theta> bv_r_min(r_min);
    ddc::ConstantExtrapolationRule<R, Theta> bv_r_max(r_max);
    ddc::PeriodicExtrapolationRule<Theta> bv_knots_min;
    ddc::PeriodicExtrapolationRule<Theta> bv_knots_max;
    SplineRThetaEvaluatorConstBound_host
            spline_evaluator_extrapol(bv_r_min, bv_r_max, bv_knots_min, bv_knots_max);

    DiscreteToCartesianBuilder<
            X,
            Y,
            SplineRThetaBuilder_host,
            SplineRThetaEvaluatorConstBound_host> const
            discrete_mapping_builder(
                    Kokkos::DefaultHostExecutionSpace(),
                    mapping_1,
                    builder,
                    spline_evaluator_extrapol);
    DiscreteToCartesian const discrete_mapping = discrete_mapping_builder();
    TOLs[0][0] = 5e-6;
    TOLs[0][1] = 5e-7;
    TOLs[1][0] = 5e-3;
    TOLs[1][1] = 1e-6;
    TOLs[2][0] = 1e-6;
    TOLs[2][1] = 5e-7;
    TOLs[3][0] = 1e-3;
    TOLs[3][1] = 5e-7;
    TOLs[4][0] = 5e-2;
    TOLs[4][1] = 5e-6;

    launch_tests(discrete_mapping, builder, grid, expected_norms, TOLs);
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        SplineQuadrature,
        testing::Combine(testing::Values<std::size_t>(40), testing::Values<std::size_t>(80)));
