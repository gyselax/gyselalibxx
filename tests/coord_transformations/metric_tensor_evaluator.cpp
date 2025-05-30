/// Test of the metric tensor and its inverse: (singular point avoided)
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "coord_transformations_testing_tools.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "geometry_coord_transformations_tests.hpp"
#include "metric_tensor_evaluator.hpp"
#include "toroidal_to_cylindrical.hpp"
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

TEST_P(InverseMetricTensor3D, InverseMatrixToroidalMap)
{
    auto const [Nrho, Ntheta, Nphi] = GetParam();
    using Mapping2D = CircularToCartesian<Rho, Theta, R, Z>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;
    using Mapping = CombinedMapping<CylindricalMapping, ToroidalMapping, Coord<Rho, Theta, Phi>>;
    double major_radius = 6.2;
    Mapping2D polar_to_RZ(major_radius);
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    Mapping mapping(cylindrical_to_cartesian, toroidal_to_cylindrical);

    FieldMemRhoThetaPhi_host<Coord<Rho, Theta, Phi>> coords
            = get_example_coords(IdxStepRho(Nrho), IdxStepTheta(Ntheta), IdxStepPhi(Nphi));
    // Exclude the centre point where the inversion is singular
    IdxStepRhoThetaPhi skip_index_step(IdxStepRho(1), IdxStepTheta(0), IdxStepPhi(0));
    IdxRangeRhoThetaPhi grid = get_idx_range(coords).remove_first(skip_index_step);

    MetricTensorEvaluator<Mapping, Coord<Rho, Theta, Phi>> metric_tensor(mapping);
    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRhoThetaPhi const idx) {
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
