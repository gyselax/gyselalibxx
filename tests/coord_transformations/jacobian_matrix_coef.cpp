#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "coord_transformations_testing_tools.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "geometry_coord_transformations_tests.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "linear_coord_transform.hpp"
#include "mesh_builder.hpp"
#include "orthogonal_coord_transforms.hpp"



/**
 * @brief A class for the Google tests.
 */
class JacobianMatrixAndJacobianCoefficients
    : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};


// Circular mapping ------------------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_jacobian_v<decltype(mapping)>);
    static_assert(has_inv_jacobian_v<decltype(mapping)>);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        Tensor Jacobian_matrix = mapping.jacobian_matrix(coords(irtheta));

        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> Jacobian_matrix_coeff;
        ddcHelper::get<X, R_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<X, R_cov>(coords(irtheta));
        ddcHelper::get<X, Theta_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<X, Theta_cov>(coords(irtheta));
        ddcHelper::get<Y, R_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<Y, R_cov>(coords(irtheta));
        ddcHelper::get<Y, Theta_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<Y, Theta_cov>(coords(irtheta));

        EXPECT_NEAR(
                (ddcHelper::get<X, R_cov>(Jacobian_matrix)),
                (ddcHelper::get<X, R_cov>(Jacobian_matrix_coeff)),
                1e-13);
        EXPECT_NEAR(
                (ddcHelper::get<X, Theta_cov>(Jacobian_matrix)),
                (ddcHelper::get<X, Theta_cov>(Jacobian_matrix_coeff)),
                1e-13);
        EXPECT_NEAR(
                (ddcHelper::get<Y, R_cov>(Jacobian_matrix)),
                (ddcHelper::get<Y, R_cov>(Jacobian_matrix_coeff)),
                1e-13);
        EXPECT_NEAR(
                (ddcHelper::get<Y, Theta_cov>(Jacobian_matrix)),
                (ddcHelper::get<Y, Theta_cov>(Jacobian_matrix_coeff)),
                1e-13);
    });

    // --- for the inverse Jacobian matrix:
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        Tensor inv_Jacobian_matrix = mapping.inv_jacobian_matrix(coords(irtheta));

        DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X, Y>> inv_Jacobian_matrix_coeff;
        ddcHelper::get<R, X>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<R, X>(coords(irtheta));
        ddcHelper::get<R, Y>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<R, Y>(coords(irtheta));
        ddcHelper::get<Theta, X>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<Theta, X>(coords(irtheta));
        ddcHelper::get<Theta, Y>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<Theta, Y>(coords(irtheta));

        double mat_elem = ddcHelper::get<R, X>(inv_Jacobian_matrix);
        double coeff = ddcHelper::get<R, X>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
        mat_elem = ddcHelper::get<R, Y>(inv_Jacobian_matrix);
        coeff = ddcHelper::get<R, Y>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
        mat_elem = ddcHelper::get<Theta, X>(inv_Jacobian_matrix);
        coeff = ddcHelper::get<Theta, X>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
        mat_elem = ddcHelper::get<Theta, Y>(inv_Jacobian_matrix);
        coeff = ddcHelper::get<Theta, Y>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
    });
}



// Czarny mapping --------------------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 1.4);

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_jacobian_v<decltype(mapping)>);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        Tensor Jacobian_matrix = mapping.jacobian_matrix(coords(irtheta));

        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> Jacobian_matrix_coeff;
        ddcHelper::get<X, R_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<X, R_cov>(coords(irtheta));
        ddcHelper::get<X, Theta_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<X, Theta_cov>(coords(irtheta));
        ddcHelper::get<Y, R_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<Y, R_cov>(coords(irtheta));
        ddcHelper::get<Y, Theta_cov>(Jacobian_matrix_coeff)
                = mapping.template jacobian_component<Y, Theta_cov>(coords(irtheta));

        EXPECT_TRUE(Jacobian_matrix == Jacobian_matrix_coeff);
        EXPECT_NEAR(determinant(Jacobian_matrix), mapping.jacobian(coords(irtheta)), 1e-14);
    });

    InverseJacobianMatrix inv_jacobian(mapping);
    // --- for the inverseJacobian matrix:
    using X_cov = X::Dual;
    using Y_cov = Y::Dual;
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        Tensor inv_Jacobian_matrix = inv_jacobian(coords(irtheta));

        DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>> inv_Jacobian_matrix_coeff;
        ddcHelper::get<R, X_cov>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<R, X_cov>(coords(irtheta));
        ddcHelper::get<R, Y_cov>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<R, Y_cov>(coords(irtheta));
        ddcHelper::get<Theta, X_cov>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<Theta, X_cov>(coords(irtheta));
        ddcHelper::get<Theta, Y_cov>(inv_Jacobian_matrix_coeff)
                = mapping.template inv_jacobian_component<Theta, Y_cov>(coords(irtheta));

        double mat_elem = ddcHelper::get<R, X>(inv_Jacobian_matrix);
        double coeff = ddcHelper::get<R, X>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
        mat_elem = ddcHelper::get<R, Y>(inv_Jacobian_matrix);
        coeff = ddcHelper::get<R, Y>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
        mat_elem = ddcHelper::get<Theta, X>(inv_Jacobian_matrix);
        coeff = ddcHelper::get<Theta, X>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
        mat_elem = ddcHelper::get<Theta, Y>(inv_Jacobian_matrix);
        coeff = ddcHelper::get<Theta, Y>(inv_Jacobian_matrix_coeff);
        EXPECT_DOUBLE_EQ(mat_elem, coeff);
    });
}



// Discrete Czarny mapping -----------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixDiscCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> analytical_mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    std::vector<CoordR> r_break_points = build_uniform_break_points(r_min, r_max, r_size);
    std::vector<CoordTheta> theta_break_points
            = build_uniform_break_points(theta_min, theta_max, theta_size);

    ddc::init_discrete_space<BSplinesR>(r_break_points);
    ddc::init_discrete_space<BSplinesTheta>(theta_break_points);

    ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(InterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR interpolation_idx_range_r(InterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_theta(InterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_r, interpolation_idx_range_theta);

    SplineRThetaBuilder_host builder(grid);
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluator_host evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);
    DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator_host>
            mapping_builder(
                    Kokkos::DefaultHostExecutionSpace(),
                    analytical_mapping,
                    builder,
                    evaluator);
    DiscreteToCartesian mapping = mapping_builder();

    static_assert(has_jacobian_v<decltype(mapping)>);
    InverseJacobianMatrix inv_jacobian(mapping);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    ddc::host_for_each(grid, [&](IdxRTheta const irtheta) {
        const CoordRTheta coord_rtheta(ddc::coordinate(irtheta));
        const double r = ddc::get<R>(coord_rtheta);
        if (fabs(r) > 1e-15) {
            // --- for the Jacobian matrix:
            Tensor Jacobian_matrix = mapping.jacobian_matrix(coord_rtheta);

            DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> Jacobian_matrix_coeff;
            ddcHelper::get<X, R_cov>(Jacobian_matrix_coeff)
                    = mapping.template jacobian_component<X, R_cov>(coord_rtheta);
            ddcHelper::get<X, Theta_cov>(Jacobian_matrix_coeff)
                    = mapping.template jacobian_component<X, Theta_cov>(coord_rtheta);
            ddcHelper::get<Y, R_cov>(Jacobian_matrix_coeff)
                    = mapping.template jacobian_component<Y, R_cov>(coord_rtheta);
            ddcHelper::get<Y, Theta_cov>(Jacobian_matrix_coeff)
                    = mapping.template jacobian_component<Y, Theta_cov>(coord_rtheta);

            EXPECT_NEAR(
                    (ddcHelper::get<X, R_cov>(Jacobian_matrix)),
                    (ddcHelper::get<X, R_cov>(Jacobian_matrix_coeff)),
                    1e-13);
            EXPECT_NEAR(
                    (ddcHelper::get<X, Theta_cov>(Jacobian_matrix)),
                    (ddcHelper::get<X, Theta_cov>(Jacobian_matrix_coeff)),
                    1e-13);
            EXPECT_NEAR(
                    (ddcHelper::get<Y, R_cov>(Jacobian_matrix)),
                    (ddcHelper::get<Y, R_cov>(Jacobian_matrix_coeff)),
                    1e-13);
            EXPECT_NEAR(
                    (ddcHelper::get<Y, Theta_cov>(Jacobian_matrix)),
                    (ddcHelper::get<Y, Theta_cov>(Jacobian_matrix_coeff)),
                    1e-13);

            // --- for the inverse Jacobian matrix:
            Tensor inv_Jacobian_matrix = inv_jacobian(coord_rtheta);

            using X_cov = X::Dual;
            using Y_cov = Y::Dual;
            DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<X_cov, Y_cov>>
                    inv_Jacobian_matrix_coeff;
            ddcHelper::get<R, X_cov>(inv_Jacobian_matrix_coeff)
                    = inv_jacobian.template inv_jacobian_component<R, X_cov>(coord_rtheta);
            ddcHelper::get<R, Y_cov>(inv_Jacobian_matrix_coeff)
                    = inv_jacobian.template inv_jacobian_component<R, Y_cov>(coord_rtheta);
            ddcHelper::get<Theta, X_cov>(inv_Jacobian_matrix_coeff)
                    = inv_jacobian.template inv_jacobian_component<Theta, X_cov>(coord_rtheta);
            ddcHelper::get<Theta, Y_cov>(inv_Jacobian_matrix_coeff)
                    = inv_jacobian.template inv_jacobian_component<Theta, Y_cov>(coord_rtheta);

            EXPECT_NEAR(
                    (ddcHelper::get<R, X_cov>(inv_Jacobian_matrix)),
                    (ddcHelper::get<R, X_cov>(inv_Jacobian_matrix_coeff)),
                    1e-13);
            EXPECT_NEAR(
                    (ddcHelper::get<R, Y_cov>(inv_Jacobian_matrix)),
                    (ddcHelper::get<R, Y_cov>(inv_Jacobian_matrix_coeff)),
                    1e-13);
            EXPECT_NEAR(
                    (ddcHelper::get<Theta, X_cov>(inv_Jacobian_matrix)),
                    (ddcHelper::get<Theta, X_cov>(inv_Jacobian_matrix_coeff)),
                    1e-13);
            EXPECT_NEAR(
                    (ddcHelper::get<Theta, Y_cov>(inv_Jacobian_matrix)),
                    (ddcHelper::get<Theta, Y_cov>(inv_Jacobian_matrix_coeff)),
                    1e-13);
        }
    });
}

// Discrete Czarny mapping -----------------------------------------
TEST(JacobianMatrixAndJacobianCoefficients, OrthogonalCoordTransforms)
{
    double x_scaling(-0.5);
    double y_scaling(3.5);
    const LinearCoordTransform<X, X2> x_transform(Coord<X>(0.5), Coord<X2>(22.0), x_scaling);
    const LinearCoordTransform<Y, Y2> y_transform(Coord<Y>(-7), Coord<Y2>(0.0), y_scaling);
    OrthogonalCoordTransforms<
            Coord<X, Y>,
            Coord<Y2, X2>,
            LinearCoordTransform<X, X2>,
            LinearCoordTransform<Y, Y2>>
            transform(x_transform, y_transform);

    Coord<X, Y> test_coord(0.0, 0.0);
    Coord<Y2, X2> transformed_coord = transform(test_coord);
    EXPECT_DOUBLE_EQ(Coord<Y2>(transformed_coord), y_transform(Coord<Y>(test_coord)));
    EXPECT_DOUBLE_EQ(Coord<Y2>(transformed_coord), transform(Coord<Y>(test_coord)));
    EXPECT_DOUBLE_EQ(Coord<X2>(transformed_coord), x_transform(Coord<X>(test_coord)));
    EXPECT_DOUBLE_EQ(Coord<X2>(transformed_coord), transform(Coord<X>(test_coord)));

    EXPECT_DOUBLE_EQ((transform.template jacobian_component<X2, X>(test_coord)), x_scaling);
    EXPECT_DOUBLE_EQ((transform.template jacobian_component<Y2, Y>(test_coord)), y_scaling);
    EXPECT_DOUBLE_EQ((transform.template jacobian_component<X2, Y>(test_coord)), 0.0);
    EXPECT_DOUBLE_EQ((transform.template jacobian_component<Y2, X>(test_coord)), 0.0);

    DTensor<VectorIndexSet<Y2, X2>, VectorIndexSet<X, Y>> Jacobian_matrix
            = transform.jacobian_matrix(test_coord);
    EXPECT_DOUBLE_EQ((ddcHelper::get<X2, X>(Jacobian_matrix)), x_scaling);
    EXPECT_DOUBLE_EQ((ddcHelper::get<Y2, Y>(Jacobian_matrix)), y_scaling);
    EXPECT_DOUBLE_EQ((ddcHelper::get<X2, Y>(Jacobian_matrix)), 0.0);
    EXPECT_DOUBLE_EQ((ddcHelper::get<Y2, X>(Jacobian_matrix)), 0.0);
    EXPECT_DOUBLE_EQ(transform.jacobian(test_coord), -x_scaling * y_scaling);
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        JacobianMatrixAndJacobianCoefficients,
        testing::Combine(testing::Values<std::size_t>(40), testing::Values<std::size_t>(80)));
