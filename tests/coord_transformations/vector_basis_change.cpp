// SPDX-License-Identifier: MIT
#include <vector>

#include <gtest/gtest.h>

#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "geometry_coord_transformations_tests.hpp"
#include "mesh_builder.hpp"
#include "species_info.hpp"
#include "vector_mapper.hpp"

namespace {
static constexpr int N_cells_r = 10;
static constexpr int N_cells_theta = 10;
void setup_coordinates()
{
    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_ncells(N_cells_r);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_ncells(N_cells_theta);

    std::vector<CoordR> r_break_points = build_uniform_break_points(r_min, r_max, r_ncells);
    std::vector<CoordTheta> theta_break_points
            = build_uniform_break_points(theta_min, theta_max, theta_ncells);

    ddc::init_discrete_space<BSplinesR>(r_break_points);
    ddc::init_discrete_space<BSplinesTheta>(theta_break_points);

    ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(InterpPointsTheta::get_sampling<GridTheta>());
}
}; // namespace

TEST(MappingChange, ArgToResult)
{
    setup_coordinates();
    const CartesianToCircular<X, Y, R, Theta> mapping;
    DVector<X, Y> A_cart(1., 1.);
    Coord<X, Y> test_coord(1.0, 1.0);
    DVector<R, Theta> A_polar
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, test_coord, A_cart);
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
    Coord<R, Theta> test_coord(std::sqrt(2), 0.25 * M_PI);
    DVector<R, Theta> A_polar
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, test_coord, A_cart);
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
    Coord<R, Theta> test_coord(std::sqrt(2), 0.25 * M_PI);
    DVector<R_cov, Theta_cov> A_cov
            = to_vector_space<VectorIndexSet<R_cov, Theta_cov>>(mapping, test_coord, A_contra);
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
    Coord<R, Theta> test_coord(std::sqrt(2), 0.25 * M_PI);
    DVector<R, Theta> A_contra
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, test_coord, A_cov);
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
    Coord<R, Theta> test_coord(std::sqrt(2), 0.25 * M_PI);
    DVector<X, Y> A_cart = to_vector_space<VectorIndexSet<X, Y>>(mapping, test_coord, A_polar_cov);
    double const x = ddcHelper::get<X>(A_cart);
    double const y = ddcHelper::get<Y>(A_cart);
    DVector<R, Theta> A_polar_contra
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, test_coord, A_polar_cov);
    DVector<X, Y> A_cart_from_contra
            = to_vector_space<VectorIndexSet<X, Y>>(mapping, test_coord, A_polar_contra);
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
    Coord<R, Theta> test_coord(std::sqrt(2), 0.25 * M_PI);
    DVector<R_cov, Theta_cov> A_polar_cov
            = to_vector_space<VectorIndexSet<R_cov, Theta_cov>>(mapping, test_coord, A_cart);
    double const r_cov = ddcHelper::get<R_cov>(A_polar_cov);
    double const theta_cov = ddcHelper::get<Theta_cov>(A_polar_cov);
    DVector<R, Theta> A_polar_contra
            = to_vector_space<VectorIndexSet<R, Theta>>(mapping, test_coord, A_cart);
    DVector<R_cov, Theta_cov> A_polar_cov_from_contra = to_vector_space<
            VectorIndexSet<R_cov, Theta_cov>>(mapping, test_coord, A_polar_contra);
    double const r_cov_from_contra = ddcHelper::get<R_cov>(A_polar_cov_from_contra);
    double const theta_cov_from_contra = ddcHelper::get<Theta_cov>(A_polar_cov_from_contra);
    EXPECT_NEAR(r_cov, r_cov_from_contra, 1e-14);
    EXPECT_NEAR(theta_cov, theta_cov_from_contra, 1e-14);
}

TEST(MappingChange, VectorFieldChange)
{
    setup_coordinates();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    IdxRange<Species> idx_range_sp(Idx<Species>(0), IdxStep<Species>(3));
    IdxRangeR idx_range_r(IdxR(1), IdxStepR(N_cells_r)); // Exclude r = 0
    IdxRangeTheta idx_range_theta(IdxTheta(0), IdxStepTheta(N_cells_theta));
    IdxRange<Species, GridR, GridTheta> idx_range(idx_range_sp, idx_range_r, idx_range_theta);

    using CartesianBasis = VectorIndexSet<X, Y>;
    using PolarBasis = VectorIndexSet<R, Theta>;

    DVectorFieldMem<IdxRange<Species, GridR, GridTheta>, CartesianBasis> cart_vector_field(
            idx_range);
    ddc::parallel_fill(ddcHelper::get<X>(cart_vector_field), 1.0);
    ddc::parallel_fill(ddcHelper::get<Y>(cart_vector_field), 1.0);

    auto polar_vector_field = create_mirror_view_and_copy_on_vector_space<
            PolarBasis>(Kokkos::DefaultExecutionSpace(), get_field(cart_vector_field), mapping);

    auto polar_vector_field_host = ddcHelper::create_mirror_view_and_copy(
            Kokkos::DefaultHostExecutionSpace(),
            get_field(polar_vector_field));

    ddc::for_each(idx_range, [&](Idx<Species, GridR, GridTheta> idx) {
        CoordRTheta coord = ddc::coordinate(IdxRTheta(idx));
        Tensor J = mapping.inv_jacobian_matrix(coord);
        double r_vec = ddcHelper::get<R>(polar_vector_field_host)(idx);
        double theta_vec = ddcHelper::get<Theta>(polar_vector_field_host)(idx);
        double expected_r_vec = ddcHelper::get<R, X>(J) + ddcHelper::get<R, Y>(J);
        double expected_theta_vec = ddcHelper::get<Theta, X>(J) + ddcHelper::get<Theta, Y>(J);
        EXPECT_NEAR(r_vec, expected_r_vec, 1e-14);
        EXPECT_NEAR(theta_vec, expected_theta_vec, 1e-14);
    });
}

TEST(MappingChange, VectorFieldChangeCopy)
{
    setup_coordinates();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    IdxRange<Species> idx_range_sp(Idx<Species>(0), IdxStep<Species>(3));
    IdxRangeR idx_range_r(IdxR(1), IdxStepR(N_cells_r)); // Exclude r = 0
    IdxRangeTheta idx_range_theta(IdxTheta(0), IdxStepTheta(N_cells_theta));
    IdxRange<Species, GridR, GridTheta> idx_range(idx_range_sp, idx_range_r, idx_range_theta);

    using CartesianBasis = VectorIndexSet<X, Y>;
    using PolarBasis = VectorIndexSet<R, Theta>;

    DVectorFieldMem<IdxRange<Species, GridR, GridTheta>, CartesianBasis> cart_vector_field(
            idx_range);
    ddc::parallel_fill(ddcHelper::get<X>(cart_vector_field), 1.0);
    ddc::parallel_fill(ddcHelper::get<Y>(cart_vector_field), 1.0);

    DVectorFieldMem<IdxRange<Species, GridR, GridTheta>, PolarBasis> polar_vector_field(idx_range);
    copy_to_vector_space<PolarBasis>(
            Kokkos::DefaultExecutionSpace(),
            get_field(polar_vector_field),
            mapping,
            get_const_field(cart_vector_field));

    auto polar_vector_field_host = ddcHelper::create_mirror_view_and_copy(
            Kokkos::DefaultHostExecutionSpace(),
            get_field(polar_vector_field));

    ddc::for_each(idx_range, [&](Idx<Species, GridR, GridTheta> idx) {
        CoordRTheta coord = ddc::coordinate(IdxRTheta(idx));
        Tensor J = mapping.inv_jacobian_matrix(coord);
        double r_vec = ddcHelper::get<R>(polar_vector_field_host)(idx);
        double theta_vec = ddcHelper::get<Theta>(polar_vector_field_host)(idx);
        double expected_r_vec = ddcHelper::get<R, X>(J) + ddcHelper::get<R, Y>(J);
        double expected_theta_vec = ddcHelper::get<Theta, X>(J) + ddcHelper::get<Theta, Y>(J);
        EXPECT_NEAR(r_vec, expected_r_vec, 1e-14);
        EXPECT_NEAR(theta_vec, expected_theta_vec, 1e-14);
    });
}
