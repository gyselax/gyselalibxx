// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "../mapping/geometry_mapping_tests.hpp"

#include "circular_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "metric_tensor_evaluator.hpp"
#include "tensor.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

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
