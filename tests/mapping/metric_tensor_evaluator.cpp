/// Test of the metric tensor and its inverse: (singular point avoided)
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "geometry_mapping_tests.hpp"
#include "mapping_testing_tools.hpp"
#include "metric_tensor_evaluator.hpp"
#include "view.hpp"



class InverseMetricTensor : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};

class InverseMetricTensor3D
    : public testing::TestWithParam<std::tuple<std::size_t, std::size_t, std::size_t>>
{
};

TEST_P(InverseMetricTensor, InverseMatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    MetricTensorEvaluator<CircularToCartesian<R, Theta, X, Y>, CoordRTheta> metric_tensor(mapping);
    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        check_inverse_tensor(
                metric_tensor(coords(irtheta)),
                metric_tensor.inverse(coords(irtheta)),
                1e-10);
    });
}



TEST_P(InverseMetricTensor, InverseMatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 1.4);

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    MetricTensorEvaluator<CzarnyToCartesian<R, Theta, X, Y>, CoordRTheta> metric_tensor(mapping);
    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        check_inverse_tensor(
                metric_tensor(coords(irtheta)),
                metric_tensor.inverse(coords(irtheta)),
                1e-10);
    });
}

TEST_P(InverseMetricTensor3D, InverseMatrixCylindricalMap)
{
    auto const [Nr, Nz, Nzeta] = GetParam();
    using Mapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;
    Mapping mapping;

    FieldMemRZZeta_host<Coord<R, Z, Zeta>> coords
            = get_example_coords(IdxStepR(Nr), IdxStepZ(Nz), IdxStepZeta(Nzeta));
    IdxRangeRZZeta grid = get_idx_range(coords);

    MetricTensorEvaluator<Mapping, Coord<R, Z, Zeta>> metric_tensor(mapping);
    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRZZeta const idx) {
        check_inverse_tensor(metric_tensor(coords(idx)), metric_tensor.inverse(coords(idx)), 1e-10);
    });
}

INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InverseMetricTensor,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(64)));

INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InverseMetricTensor3D,
        testing::Combine(
                testing::Values<std::size_t>(16),
                testing::Values<std::size_t>(16),
                testing::Values<std::size_t>(8)));
