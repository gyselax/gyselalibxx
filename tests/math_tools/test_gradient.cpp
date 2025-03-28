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
#include "mesh_builder.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"

namespace {

using DVectorType = DVector<R, Theta>;
using DVectorCovType = DVector<R_cov, Theta_cov>;

using DVectorFieldMemType = DVectorFieldMem<IdxRangeRTheta, VectorIndexSet<R, Theta>>;
using DVectorFieldType = typename DVectorFieldMemType::span_type;

using DVectorFieldMemCovType = DVectorFieldMem<IdxRangeRTheta, VectorIndexSet<R_cov, Theta_cov>>;
using DVectorFieldCovType = typename DVectorFieldMemCovType::span_type;

template <class ElementType>
using FieldMemRTheta = FieldMem<ElementType, IdxRangeRTheta>;

template <class ElementType>
using FieldRTheta = Field<ElementType, IdxRangeRTheta>;


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
    KOKKOS_INLINE_FUNCTION double operator()(CoordRTheta const coord_rth) const
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
    KOKKOS_FUNCTION DVectorType gradient(CoordRTheta const coord_rth) const
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
    KOKKOS_FUNCTION DVectorCovType partial_derivatives(CoordRTheta const coord_rth) const
    {
        double const r = ddc::get<R>(coord_rth);
        double const theta = ddc::get<Theta>(coord_rth);
        DVectorCovType derivatives(
                2 * r * Kokkos::cos(theta) + 6 * ipow(r, 2) * Kokkos::sin(theta),
                -ipow(r, 2) * Kokkos::sin(theta) + 2 * ipow(r, 3) * Kokkos::cos(theta));

        return derivatives;
    }
};

// Function to compare the gradient to its prediction elementwise
void test_gradient_prediction(
        DVectorFieldType gradient_prediction_contra,
        DVectorFieldCovType gradient_prediction_cov,
        DVectorFieldType gradient_contra,
        DVectorFieldCovType gradient_cov,
        double const TOL)
{
    auto gradient_prediction_contra_host = ddcHelper::create_mirror_view_and_copy(
            Kokkos::DefaultHostExecutionSpace(),
            gradient_prediction_contra);
    auto gradient_prediction_cov_host = ddcHelper::create_mirror_view_and_copy(
            Kokkos::DefaultHostExecutionSpace(),
            gradient_prediction_cov);
    auto gradient_contra_host = ddcHelper::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), gradient_contra);
    auto gradient_cov_host = ddcHelper::
            create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), gradient_cov);

    ddc::for_each(get_idx_range(gradient_contra_host), [&](IdxRTheta const irtheta) {
        EXPECT_NEAR(
                ddcHelper::get<R>(gradient_contra_host)(irtheta),
                ddcHelper::get<R>(gradient_prediction_contra_host)(irtheta),
                TOL);
        EXPECT_NEAR(
                ddcHelper::get<Theta>(gradient_contra_host)(irtheta),
                ddcHelper::get<Theta>(gradient_prediction_contra_host)(irtheta),
                TOL);

        EXPECT_NEAR(
                ddcHelper::get<R_cov>(gradient_cov_host)(irtheta),
                ddcHelper::get<R_cov>(gradient_prediction_cov_host)(irtheta),
                TOL);
        EXPECT_NEAR(
                ddcHelper::get<Theta_cov>(gradient_cov_host)(irtheta),
                ddcHelper::get<Theta_cov>(gradient_prediction_cov_host)(irtheta),
                TOL);
    });
}

void compute_and_test_gradient()
{
    double const TOL = 1e-12;
    CircularToCartesian<R, Theta, X, Y> const mapping;

    IdxStepR nbcells_r(32);
    std::vector<CoordR> point_sampling_r
            = build_uniform_break_points(CoordR(1e-5), CoordR(1.), nbcells_r);

    ddc::init_discrete_space<GridR>(point_sampling_r);

    IdxStepTheta nbcells_th(32);
    std::vector<CoordTheta> point_sampling_th
            = build_uniform_break_points(CoordTheta(0.), CoordTheta(2. * M_PI), nbcells_th);
    ddc::init_discrete_space<GridTheta>(point_sampling_th);

    IdxRangeR idxrange_r(IdxR(0), nbcells_r + 1);
    IdxRangeTheta idxrange_th(IdxTheta(0), nbcells_th);
    IdxRangeRTheta idxrange_rtheta(idxrange_r, idxrange_th);

    using MetricTensorType
            = MetricTensorEvaluator<CircularToCartesian<R, Theta, X, Y>, CoordRTheta>;
    MetricTensorType metric_tensor(mapping);
    Gradient<MetricTensorType> gradient_evaluator(metric_tensor);

    GradientTestFunction test_function;

    // computation of partial derivatives, gradient, and predictions for the gradient
    DVectorFieldMemCovType partial_derivatives_alloc(idxrange_rtheta);
    DVectorFieldCovType partial_derivatives = get_field(partial_derivatives_alloc);

    DVectorFieldMemType gradient_prediction_contra_alloc(idxrange_rtheta);
    DVectorFieldType gradient_prediction_contra = get_field(gradient_prediction_contra_alloc);

    DVectorFieldMemType gradient_contra_alloc(idxrange_rtheta);
    DVectorFieldType gradient_contra = get_field(gradient_contra_alloc);

    DVectorFieldMemCovType gradient_cov_alloc(idxrange_rtheta);
    DVectorFieldCovType gradient_cov = get_field(gradient_cov_alloc);

    // Test for each coordinate if the gradient is equal to its prediction
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_rtheta,
            KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                CoordRTheta const coord_rth(ddc::coordinate(irtheta));

                ddcHelper::assign_vector_field_element(
                        partial_derivatives,
                        irtheta,
                        test_function.partial_derivatives(coord_rth));

                ddcHelper::assign_vector_field_element(
                        gradient_prediction_contra,
                        irtheta,
                        test_function.gradient(coord_rth));

                ddcHelper::assign_vector_field_element(
                        gradient_contra,
                        irtheta,
                        gradient_evaluator(partial_derivatives(irtheta), coord_rth));

                ddcHelper::assign_vector_field_element(
                        gradient_cov,
                        irtheta,
                        gradient_evaluator(partial_derivatives(irtheta)));
            });

    test_gradient_prediction(
            gradient_prediction_contra,
            partial_derivatives,
            gradient_contra,
            gradient_cov,
            TOL);

    // Test for each coordinate if the gradient is equal to its prediction
    //using gradient methods that take vector field as parameters
    gradient_evaluator(gradient_contra, get_const_field(partial_derivatives));
    gradient_evaluator(gradient_cov, get_const_field(partial_derivatives));

    test_gradient_prediction(
            gradient_prediction_contra,
            partial_derivatives,
            gradient_contra,
            gradient_cov,
            TOL);
}

TEST(GradientTest, Circular)
{
    compute_and_test_gradient();
}
} // namespace
