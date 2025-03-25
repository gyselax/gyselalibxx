/// Test of the metric tensor and its inverse: (singular point avoided)
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "../mapping/geometry_mapping_tests.hpp"
#include "../mapping/mapping_testing_tools.hpp"

#include "circular_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "gradient.hpp"
#include "math_tools.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"

namespace {

using DVectorType = DTensor<VectorIndexSet<R, Theta>>;
using DVectorCovType = DTensor<VectorIndexSet<R_cov, Theta_cov>>;


/**
 * @brief A class that represents a scalar field.
 * To be used to compute the gradient of the field.
 */
class GradientTestFunction
{
public:
    /**
     * @brief Default constructor.
     */
    KOKKOS_DEFAULTED_FUNCTION GradientTestFunction() = default;

    /**
     * @brief Default copy constructor.
     */
    KOKKOS_DEFAULTED_FUNCTION GradientTestFunction(GradientTestFunction const&) = default;

    /**
     * @brief Default destructor.
     */
    KOKKOS_DEFAULTED_FUNCTION ~GradientTestFunction() = default;

    /**
     * @brief Get the value of the function at given coordinate.
     *
     * @param[in] coord_rth The coordinate where we want to evaluate
     * the function.
     *
     * @return The value of the function at the coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordRTheta const coord_rth) const
    {
        double const r = ddc::get<R>(coord_rth);
        double const theta = ddc::get<Theta>(coord_rth);
        return ipow(r, 2) * Kokkos::cos(theta) + 2 * ipow(r, 3) * Kokkos::sin(theta);
    }

    /**
     * @brief Get the value of the gradient of the field 
     * at a given coordinate. The gradient is expressed in 
     * the contravariant basis of the circular geometry. 
     *
     * @param[in] coord_rth The coordinate where we want to evaluate 
     * the gradient.
     *
     * @return A vector type representing the value of the gradient.
     */
    KOKKOS_INLINE_FUNCTION DVectorType gradient(CoordRTheta const coord_rth) const
    {
        double const r = ddc::get<R>(coord_rth);
        double const theta = ddc::get<Theta>(coord_rth);
        DVectorType gradient(
                2 * r * Kokkos::cos(theta) + 6 * ipow(r, 2) * Kokkos::sin(theta),
                -Kokkos::sin(theta) + 2 * r * Kokkos::cos(theta));

        return gradient;
    }

    /**
     * @brief Get the value of the partial derivatives of the field 
     * at a given coordinate.
     *
     * @param[in] coord_rth The coordinate where we want to evaluate 
     * the partial derivatives.
     *
     * @return A vector type representing the value of the partial derivatives.
     */
    KOKKOS_INLINE_FUNCTION DVectorCovType partial_derivatives(CoordRTheta const coord_rth) const
    {
        double const r = ddc::get<R>(coord_rth);
        double const theta = ddc::get<Theta>(coord_rth);
        DVectorCovType derivatives(
                2 * r * Kokkos::cos(theta) + 6 * ipow(r, 2) * Kokkos::sin(theta),
                -ipow(r, 2) * Kokkos::sin(theta) + 2 * ipow(r, 3) * Kokkos::cos(theta));
        return derivatives;
    }
};


TEST(GradientTest, Circular)
{
    double const TOL = 1e-12;
    CircularToCartesian<R, Theta, X, Y> const mapping;
    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(32), IdxStepTheta(32));
    IdxRangeRTheta grid = get_idx_range(coords);

    MetricTensorEvaluator<CircularToCartesian<R, Theta, X, Y>, CoordRTheta> metric_tensor(mapping);
    Gradient<CircularToCartesian<R, Theta, X, Y>, CoordRTheta> gradient_evaluator(metric_tensor);

    GradientTestFunction test_function;

    // Test for each coordinates if the gradient is equal to its prediction
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        CoordRTheta const coord_rth(coords(irtheta));
        DVectorCovType partial_derivatives = test_function.partial_derivatives(coord_rth);

        // contravariant expression of the gradient
        DVectorType gradient_contra = gradient_evaluator(partial_derivatives, coord_rth);
        DVectorType gradient_prediction_contra = test_function.gradient(coord_rth);

        EXPECT_NEAR(
                ddcHelper::get<R>(gradient_contra),
                ddcHelper::get<R>(gradient_prediction_contra),
                TOL);
        EXPECT_NEAR(
                ddcHelper::get<Theta>(gradient_contra),
                ddcHelper::get<Theta>(gradient_prediction_contra),
                TOL);

        // covariant expression of the gradient
        DVectorCovType gradient_cov = gradient_evaluator(partial_derivatives);

        EXPECT_NEAR(
                ddcHelper::get<R_cov>(gradient_cov),
                ddcHelper::get<R_cov>(partial_derivatives),
                TOL);
        EXPECT_NEAR(
                ddcHelper::get<Theta_cov>(gradient_cov),
                ddcHelper::get<Theta_cov>(partial_derivatives),
                TOL);
    });
}
} // namespace
