// SPDX-License-Identifier: MIT
#include <vector>

#include <gtest/gtest.h>

#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "geometry_mapping_tests.hpp"
#include "mesh_builder.hpp"
#include "vector_mapper.hpp"

namespace {
void setup_coordinates()
{
    int Nr = 10;
    int Ntheta = 10;
    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Ntheta);

    std::vector<CoordR> r_knots = build_uniform_break_points(r_min, r_max, r_size);
    std::vector<CoordTheta> theta_knots
            = build_uniform_break_points(theta_min, theta_max, theta_size);

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(InterpPointsTheta::get_sampling<GridTheta>());
}
}; // namespace

TEST(MappingChange, ArgToResult)
{
    setup_coordinates();
    const CartesianToCircular<X, Y, R, Theta> mapping;
    DVector<X, Y> A_cart(1., 1.);
    Coord<X, Y> coord_centre(1.0, 1.0);
    DVector<R, Theta> A_polar
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, coord_centre, A_cart);
    double const r = ddcHelper::get<R>(A_polar);
    double const theta = ddcHelper::get<Theta>(A_polar);
    EXPECT_NEAR(r, std::sqrt(2), 1e-14);
    EXPECT_NEAR(theta, 0., 1e-14);
}

TEST(MappingChange, ResultToArg)
{
    setup_coordinates();
    const CircularToCartesian<R, Theta, X, Y> mapping;
    DVector<X, Y> A_cart(1., 1.);
    Coord<R, Theta> coord_centre(std::sqrt(2), 0.25 * M_PI);
    DVector<R, Theta> A_polar
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, coord_centre, A_cart);
    double const r = ddcHelper::get<R>(A_polar);
    double const theta = ddcHelper::get<Theta>(A_polar);
    EXPECT_NEAR(r, std::sqrt(2), 1e-14);
    EXPECT_NEAR(theta, 0., 1e-14);
}

TEST(MappingChange, ContraToCov)
{
    setup_coordinates();
    const CircularToCartesian<R, Theta, X, Y> mapping;
    DVector<R, Theta> A_contra(1., 1.);
    Coord<R, Theta> coord_centre(std::sqrt(2), 0.25 * M_PI);
    DVector<R_cov, Theta_cov> A_cov
            = to_vector_space<VectorIndexSet<R_cov, Theta_cov>>(mapping, coord_centre, A_contra);
    double const r_cov = ddcHelper::get<R_cov>(A_cov);
    double const theta_cov = ddcHelper::get<Theta_cov>(A_cov);
    EXPECT_NEAR(r_cov, 1., 1e-14);
    EXPECT_NEAR(theta_cov, 2., 1e-14);
}

TEST(MappingChange, CovToContra)
{
    setup_coordinates();
    const CircularToCartesian<R, Theta, X, Y> mapping;
    DVector<R_cov, Theta_cov> A_cov(1., 1.);
    Coord<R, Theta> coord_centre(std::sqrt(2), 0.25 * M_PI);
    DVector<R, Theta> A_contra
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, coord_centre, A_cov);
    double const r = ddcHelper::get<R>(A_contra);
    double const theta = ddcHelper::get<Theta>(A_contra);
    EXPECT_NEAR(r, 1., 1e-14);
    EXPECT_NEAR(theta, 0.5, 1e-14);
}

TEST(MappingChange, CovArgToCovResult)
{
    setup_coordinates();
    const CircularToCartesian<R, Theta, X, Y> mapping;
    DVector<R_cov, Theta_cov> A_polar_cov(std::sqrt(8), 0.0);
    Coord<R, Theta> coord_centre(std::sqrt(2), 0.25 * M_PI);
    DVector<X, Y> A_cart
            = to_vector_space<VectorIndexSet<X, Y>>(mapping, coord_centre, A_polar_cov);
    double const x = ddcHelper::get<X>(A_cart);
    double const y = ddcHelper::get<Y>(A_cart);
    DVector<R, Theta> A_polar_contra
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, coord_centre, A_polar_cov);
    DVector<X, Y> A_cart_from_contra
            = to_vector_space<VectorIndexSet<X, Y>>(mapping, coord_centre, A_polar_contra);
    double const x_from_contra = ddcHelper::get<X>(A_cart_from_contra);
    double const y_from_contra = ddcHelper::get<Y>(A_cart_from_contra);
    EXPECT_NEAR(x, x_from_contra, 1e-14);
    EXPECT_NEAR(y, y_from_contra, 1e-14);
}

TEST(MappingChange, CovResultToCovArg)
{
    setup_coordinates();
    const CircularToCartesian<R, Theta, X, Y> mapping;
    DVector<X, Y> A_cart(1., 1.);
    Coord<R, Theta> coord_centre(std::sqrt(2), 0.25 * M_PI);
    DVector<R_cov, Theta_cov> A_polar_cov
            = to_vector_space<VectorIndexSet<R_cov, Theta_cov>>(mapping, coord_centre, A_cart);
    double const r_cov = ddcHelper::get<R_cov>(A_polar_cov);
    double const theta_cov = ddcHelper::get<Theta_cov>(A_polar_cov);
    DVector<R, Theta> A_polar_contra
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, coord_centre, A_cart);
    DVector<R_cov, Theta_cov> A_polar_cov_from_contra = to_vector_space<
            VectorIndexSet<R_cov, Theta_cov>>(mapping, coord_centre, A_polar_contra);
    double const r_cov_from_contra = ddcHelper::get<R_cov>(A_polar_cov_from_contra);
    double const theta_cov_from_contra = ddcHelper::get<Theta_cov>(A_polar_cov_from_contra);
    EXPECT_NEAR(r_cov, r_cov_from_contra, 1e-14);
    EXPECT_NEAR(theta_cov, theta_cov_from_contra, 1e-14);
}
