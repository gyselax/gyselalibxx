// SPDX-License-Identifier: MIT
//
// Tests the GMGPolarPoissonLikeSolver on a circular domain using a manufactured solution.
//
// Exact solution:  phi(r, theta) = C * r^6 * (r-1)^6 * cos(3*theta), C = 1e-4 * 2^12
//
// With alpha=1, beta=0, the circular mapping gives the standard polar Laplacian:
//   -Delta phi = -(phi_rr + phi_r/r + phi_tt/r^2) = rho
//
// The test passes if the L-inf error is below 1e-2.

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "geometry_r_theta.hpp"
#include "gmg_polar_poisson_like_solver.hpp"
#include "mesh_builder.hpp"
#include "spline_definitions_r_theta.hpp"

using Mapping = CircularToCartesian<R, Theta, X, Y>;

using GMGSolver = GMGPolarPoissonLikeSolver<
        Mapping,
        GridR,
        GridTheta,
        BSplinesR,
        BSplinesTheta,
        SplineRThetaBuilder_host,
        SplineRThetaEvaluatorNullBound_host>;

namespace {

// C = 1e-4 * 2^12
constexpr double C_MANUF = 1e-4 * 4096.0;
// Angular mode number
constexpr int M_MODE = 3;

double phi_exact(double r, double theta)
{
    return C_MANUF * std::pow(r, 6) * std::pow(r - 1.0, 6) * std::cos(M_MODE * theta);
}

// Analytical RHS: rho = -(f'' + f'/r - m^2 * f / r^2) * cos(m*theta)
// where f(r) = C * r^6 * (r-1)^6
double rho_exact(double r, double theta)
{
    if (r == 0.0) {
        return 0.0;
    }
    const double f = C_MANUF * std::pow(r, 6) * std::pow(r - 1.0, 6);
    const double fp = 6.0 * C_MANUF * std::pow(r, 5) * std::pow(r - 1.0, 5) * (2.0 * r - 1.0);
    const double fpp = 6.0 * C_MANUF * std::pow(r, 4) * std::pow(r - 1.0, 4)
                       * (5.0 * std::pow(2.0 * r - 1.0, 2) + 2.0 * r * (r - 1.0));
    return -(fpp + fp / r - static_cast<double>(M_MODE * M_MODE) * f / (r * r))
           * std::cos(M_MODE * theta);
}

} // namespace

TEST(GMGPolarIntegration, PoissonEquation)
{
    CoordR const r_min(1e-5);
    CoordR const r_max(1.0);
    IdxStepR const r_ncells(16 - BSDegreeR + 1); // To have a number of grid cells that is 2**N

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(64);

    std::vector<CoordR> r_break_points = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> theta_break_points
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);

    ddc::init_discrete_space<BSplinesR>(r_break_points);
    ddc::init_discrete_space<BSplinesTheta>(theta_break_points);
    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR const idx_range_r(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta const idx_range_theta(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta const grid(idx_range_r, idx_range_theta);

    ddc::NullExtrapolationRule bv_r_min;
    ddc::NullExtrapolationRule bv_r_max;
    ddc::PeriodicExtrapolationRule<Theta> bv_theta_min;
    ddc::PeriodicExtrapolationRule<Theta> bv_theta_max;

    SplineRThetaBuilder_host const builder_host(grid);
    SplineRThetaEvaluatorNullBound_host const
            evaluator_host(bv_r_min, bv_r_max, bv_theta_min, bv_theta_max);

    const Mapping mapping;

    GMGSolver solver(mapping, builder_host, evaluator_host);

    // Set alpha = 1, beta = 0 everywhere
    host_t<DFieldMemRTheta> alpha_host(grid);
    host_t<DFieldMemRTheta> beta_host(grid);
    ddc::parallel_fill(get_field(alpha_host), 1.0);
    ddc::parallel_fill(get_field(beta_host), 0.0);

    DFieldMemRTheta alpha_alloc(grid);
    DFieldMemRTheta beta_alloc(grid);
    ddc::parallel_deepcopy(get_field(alpha_alloc), get_const_field(alpha_host));
    ddc::parallel_deepcopy(get_field(beta_alloc), get_const_field(beta_host));

    solver.update_coefficients(get_const_field(alpha_alloc), get_const_field(beta_alloc));

    // Compute the RHS on the grid
    host_t<DFieldMemRTheta> rho_host(grid);
    ddc::host_for_each(grid, [&](IdxRTheta irtheta) {
        double const r = ddc::coordinate(ddc::select<GridR>(irtheta));
        double const theta = ddc::coordinate(ddc::select<GridTheta>(irtheta));
        rho_host(irtheta) = rho_exact(r, theta);
    });

    DFieldMemRTheta rho_alloc(grid);
    ddc::parallel_deepcopy(get_field(rho_alloc), get_const_field(rho_host));

    // Solve
    DFieldMemRTheta phi_alloc(grid);
    ddc::parallel_fill(get_field(phi_alloc), 0.0);
    solver(get_field(phi_alloc), get_const_field(rho_alloc));

    // Check L-inf error
    auto phi_result_host = ddc::create_mirror_view_and_copy(get_const_field(phi_alloc));
    double max_err = 0.0;
    ddc::host_for_each(grid, [&](IdxRTheta irtheta) {
        double const r = ddc::coordinate(ddc::select<GridR>(irtheta));
        double const theta = ddc::coordinate(ddc::select<GridTheta>(irtheta));
        double const err = std::abs(phi_result_host(irtheta) - phi_exact(r, theta));
        max_err = std::max(max_err, err);
    });

    std::cout << "Max L-inf error: " << max_err << std::endl;

    EXPECT_LT(max_err, 1e-6);
}
