// SPDX-License-Identifier: MIT
#include <cmath>
#include <vector>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include "../coord_transformations/geometry_coord_transformations_tests.hpp"

#include "circular_to_cartesian.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "identity_coordinate_change.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "metric_tensor_evaluator.hpp"
#include "tensor.hpp"
#include "tensor_common.hpp"
#include "vector_field.hpp"
#include "vector_field_common.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"
#include "vector_mapper.hpp"

namespace {
void vector_field_norm_test()
{
    using IndexSet = VectorIndexSet<R, Theta>;
    double constexpr TOL(1e-14);

    CoordR r_min(0.1);
    CoordR r_max(1.0);
    CoordTheta theta_min(0.0);
    CoordTheta theta_max(6.0);
    IdxR idx_r0(0);
    IdxStepR nr(10);
    IdxTheta idx_theta0(0);
    IdxStepTheta ntheta(10);

    ddc::init_discrete_space<GridR>(
            GridR::init<GridR>(build_uniform_break_points(r_min, r_max, nr)));
    ddc::init_discrete_space<GridTheta>(
            GridTheta::init<GridTheta>(build_uniform_break_points(theta_min, theta_max, ntheta)));

    IdxRangeR idx_range_r(idx_r0, nr);
    IdxRangeTheta idx_range_theta(idx_theta0, ntheta);
    IdxRangeRTheta idx_range(idx_range_r, idx_range_theta);
    DVectorFieldMem<IdxRangeRTheta, IndexSet> vec_field_alloc(idx_range);
    DVectorField<IdxRangeRTheta, IndexSet> vec_field = get_field(vec_field_alloc);
    DFieldMem<IdxRangeRTheta> norm_vals(idx_range);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxRTheta idx) {
                ddcHelper::get<R>(vec_field)(idx) = 1.0;
                ddcHelper::get<Theta>(vec_field)(idx) = 1.5;
            });

    CircularToCartesian<R, Theta, X, Y> mapping;
    MetricTensorEvaluator get_metric(mapping);

    norm(Kokkos::DefaultExecutionSpace(),
         get_field(norm_vals),
         get_metric,
         get_const_field(vec_field));

    auto norm_vals_host = ddc::create_mirror_view_and_copy(get_field(norm_vals));

    ddc::for_each(idx_range, [&](IdxRTheta idx) {
        CoordR r = ddc::coordinate(ddc::select<GridR>(idx));
        double expected_norm = std::sqrt(1.0 + 1.5 * 1.5 * (r * r));
        ASSERT_NEAR(expected_norm, norm_vals_host(idx), TOL);
    });
}
} // namespace

TEST(MathTools, VectorNorm)
{
    double constexpr TOL(1e-14);

    // Build an arbitrary metric tensor
    using IndexSetCov = VectorIndexSet<R_cov, Theta_cov>;
    DTensor<IndexSetCov, IndexSetCov> metric_tensor;
    ddcHelper::get<R_cov, R_cov>(metric_tensor) = 1;
    ddcHelper::get<R_cov, Theta_cov>(metric_tensor) = 2;
    ddcHelper::get<Theta_cov, R_cov>(metric_tensor) = 3;
    ddcHelper::get<Theta_cov, Theta_cov>(metric_tensor) = 4;

    DVector<R, Theta> test_vec1(1.0, 0.0);
    ASSERT_NEAR(norm(metric_tensor, test_vec1), 1, TOL);

    DVector<R, Theta> test_vec2(0.0, 1.0);
    ASSERT_NEAR(norm(metric_tensor, test_vec2), 2, TOL);

    DVector<R, Theta> test_vec3(1.0, 1.0);
    ASSERT_NEAR(norm(metric_tensor, test_vec3), std::sqrt(10), TOL);
}

TEST(MathTools, VectorFieldNorm)
{
    vector_field_norm_test();
}

TEST(MathTools, Inverse)
{
    using IndexSet = VectorIndexSet<X, Y, Z>;
    DTensor<IndexSet, IndexSet> matrix;
    ddcHelper::get<X, X>(matrix) = 1;
    ddcHelper::get<X, Y>(matrix) = 2;
    ddcHelper::get<X, Z>(matrix) = 3;
    ddcHelper::get<Y, X>(matrix) = 0;
    ddcHelper::get<Y, Y>(matrix) = 1;
    ddcHelper::get<Y, Z>(matrix) = 4;
    ddcHelper::get<Z, X>(matrix) = 5;
    ddcHelper::get<Z, Y>(matrix) = 6;
    ddcHelper::get<Z, Z>(matrix) = 0;

    DTensor<IndexSet, IndexSet> inv_matrix = inverse(matrix);
    double test_val = ddcHelper::get<X, X>(inv_matrix);
    ASSERT_NEAR(test_val, -24, 1e-12);
    test_val = ddcHelper::get<X, Y>(inv_matrix);
    ASSERT_NEAR(test_val, 18, 1e-12);
    test_val = ddcHelper::get<X, Z>(inv_matrix);
    ASSERT_NEAR(test_val, 5, 1e-12);
    test_val = ddcHelper::get<Y, X>(inv_matrix);
    ASSERT_NEAR(test_val, 20, 1e-12);
    test_val = ddcHelper::get<Y, Y>(inv_matrix);
    ASSERT_NEAR(test_val, -15, 1e-12);
    test_val = ddcHelper::get<Y, Z>(inv_matrix);
    ASSERT_NEAR(test_val, -4, 1e-12);
    test_val = ddcHelper::get<Z, X>(inv_matrix);
    ASSERT_NEAR(test_val, -5, 1e-12);
    test_val = ddcHelper::get<Z, Y>(inv_matrix);
    ASSERT_NEAR(test_val, 4, 1e-12);
    test_val = ddcHelper::get<Z, Z>(inv_matrix);
    ASSERT_NEAR(test_val, 1, 1e-12);
}

TEST(MathTools, ScalarProdCart)
{
    using IndexSet = VectorIndexSet<X, Y, Z>;
    Tensor<double, IndexSet> A, B;
    ddcHelper::get<X>(A) = 3;
    ddcHelper::get<Y>(A) = 6;
    ddcHelper::get<Z>(A) = -3;
    ddcHelper::get<X>(B) = 1;
    ddcHelper::get<Y>(B) = -4;
    ddcHelper::get<Z>(B) = 2;
    EXPECT_EQ(scalar_product(A, B), -27);
    Coord<X, Y, Z> test_coord(0.0, 0.0, 0.0);
    IdentityCoordinateChange<IndexSet, IndexSet> mapping;
    MetricTensorEvaluator get_metric(mapping);
    EXPECT_DOUBLE_EQ(scalar_product(get_metric(test_coord), A, B), -27);
}

TEST(MathTools, ScalarProdCyl)
{
    using IndexSet = VectorIndexSet<R, Z, Zeta>;
    Tensor<double, IndexSet> A, B;
    ddcHelper::get<R>(A) = 3;
    ddcHelper::get<Z>(A) = 6;
    ddcHelper::get<Zeta>(A) = 6;
    ddcHelper::get<R>(B) = 1;
    ddcHelper::get<Z>(B) = 2;
    ddcHelper::get<Zeta>(B) = -4;
    Coord<R, Z, Zeta> test_coord(2.5, -4.0, 0.3);
    CylindricalToCartesian<R, Z, Zeta, X, Y> mapping;
    MetricTensorEvaluator get_metric(mapping);
    EXPECT_DOUBLE_EQ(scalar_product(get_metric(test_coord), A, B), -135);
}

TEST(MathTools, TensorProdCart)
{
    using IndexSet = VectorIndexSet<X, Y, Z>;
    Tensor<double, IndexSet> A, B, C;
    ddcHelper::get<X>(A) = 3;
    ddcHelper::get<Y>(A) = 6;
    ddcHelper::get<Z>(A) = -3;
    ddcHelper::get<X>(B) = 1;
    ddcHelper::get<Y>(B) = -4;
    ddcHelper::get<Z>(B) = 2;
    Coord<X, Y, Z> test_coord(0.0, 0.0, 0.0);
    IdentityCoordinateChange<IndexSet, IndexSet> mapping;
    C = tensor_product(mapping, test_coord, A, B);
    EXPECT_NEAR(ddcHelper::get<X>(C), 0, 1e-14);
    EXPECT_NEAR(ddcHelper::get<Y>(C), -9, 1e-14);
    EXPECT_NEAR(ddcHelper::get<Z>(C), -18, 1e-14);
}

TEST(MathTools, TensorProdCyl)
{
    using CartIndexSet = VectorIndexSet<X, Y, Z>;
    using IndexSet = VectorIndexSet<R, Z, Zeta>;

    Coord<R, Z, Zeta> test_coord(2.5, -4.0, 0.3);
    CylindricalToCartesian<R, Z, Zeta, X, Y> cyl_to_cart;
    IdentityCoordinateChange<CartIndexSet, CartIndexSet> identity_mapping;

    DTensor<IndexSet> A, B;
    DTensor<vector_index_set_dual_t<IndexSet>> C_cov;
    DTensor<CartIndexSet> A_cart, B_cart, C_cart, C_cart_via_mapping;

    ddcHelper::get<R>(A) = 3;
    ddcHelper::get<Z>(A) = 6;
    ddcHelper::get<Zeta>(A) = -3;
    ddcHelper::get<R>(B) = 1;
    ddcHelper::get<Z>(B) = -4;
    ddcHelper::get<Zeta>(B) = 2;
    A_cart = to_vector_space<CartIndexSet>(cyl_to_cart, test_coord, A);
    B_cart = to_vector_space<CartIndexSet>(cyl_to_cart, test_coord, B);

    // Calculate the tensor product in cylindrical coordinates
    C_cov = tensor_product(cyl_to_cart, test_coord, A, B);
    // Calculate the tensor product in cartesian coordinates
    C_cart = tensor_product(identity_mapping, cyl_to_cart(test_coord), A_cart, B_cart);

    // Compare the results of the calculations in different coordinate systems
    C_cart_via_mapping = to_vector_space<CartIndexSet>(cyl_to_cart, test_coord, C_cov);
    EXPECT_NEAR(ddcHelper::get<X>(C_cart_via_mapping), ddcHelper::get<X>(C_cart), 1e-14);
    EXPECT_NEAR(ddcHelper::get<Y>(C_cart_via_mapping), ddcHelper::get<Y>(C_cart), 1e-14);
    EXPECT_NEAR(ddcHelper::get<Z>(C_cart_via_mapping), ddcHelper::get<Z>(C_cart), 1e-14);
}
