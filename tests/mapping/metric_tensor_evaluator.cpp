/// Test of the metric tensor and its inverse: (singular point avoided)
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "mapping_testing_tools.hpp"
#include "metric_tensor_evaluator.hpp"
#include "view.hpp"

struct X
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = X;
};
struct Y
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Y;
};
struct R_cov;
struct Theta_cov;
struct R
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = R_cov;
};

struct Theta
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Theta_cov;
};

struct R_cov
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = R;
};

struct Theta_cov
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = Theta;
};

using CoordR = Coord<R>;
using CoordTheta = Coord<Theta>;
using CoordRTheta = Coord<R, Theta>;

int constexpr BSDegree = 3;

struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
{
};
struct BSplinesTheta : ddc::NonUniformBSplines<Theta, BSDegree>
{
};

using InterpPointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using InterpPointsTheta = ddc::GrevilleInterpolationPoints<
        BSplinesTheta,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;

struct GridR : InterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : InterpPointsTheta::interpolation_discrete_dimension_type
{
};

using IdxR = Idx<GridR>;
using IdxTheta = Idx<GridTheta>;
using IdxRTheta = Idx<GridR, GridTheta>;

using IdxStepR = IdxStep<GridR>;
using IdxStepTheta = IdxStep<GridTheta>;
using IdxStepRTheta = IdxStep<GridR, GridTheta>;

using IdxRangeRTheta = IdxRange<GridR, GridTheta>;


template <class ElementType>
using FieldMemRTheta_host = host_t<FieldMem<ElementType, IdxRangeRTheta>>;


using Matrix_2x2 = std::array<std::array<double, 2>, 2>;


namespace {

template <class Dims>
void check_inverse_tensor(
        DTensor<Dims, Dims> tensor,
        DTensor<vector_index_set_dual_t<Dims>, vector_index_set_dual_t<Dims>> inv_tensor)
{
    double TOL = 1e-10;

    using Dim0 = ddc::type_seq_element_t<0, Dims>;
    using Dim1 = ddc::type_seq_element_t<1, Dims>;
    using Dim0_cov = typename Dim0::Dual;
    using Dim1_cov = typename Dim1::Dual;

    double const id_val00
            = ddcHelper::get<Dim0, Dim0>(tensor) * ddcHelper::get<Dim0_cov, Dim0_cov>(inv_tensor)
              + ddcHelper::get<Dim0, Dim1>(tensor) * ddcHelper::get<Dim1_cov, Dim0_cov>(inv_tensor);
    EXPECT_NEAR(id_val00, 1., TOL);

    double const id_val01
            = ddcHelper::get<Dim0, Dim0>(tensor) * ddcHelper::get<Dim0_cov, Dim1_cov>(inv_tensor)
              + ddcHelper::get<Dim0, Dim1>(tensor) * ddcHelper::get<Dim1_cov, Dim1_cov>(inv_tensor);
    EXPECT_NEAR(id_val01, 0., TOL);

    double const id_val10
            = ddcHelper::get<Dim1, Dim0>(tensor) * ddcHelper::get<Dim0_cov, Dim0_cov>(inv_tensor)
              + ddcHelper::get<Dim1, Dim1>(tensor) * ddcHelper::get<Dim1_cov, Dim0_cov>(inv_tensor);
    EXPECT_NEAR(id_val10, 0., TOL);

    double const id_val11
            = ddcHelper::get<Dim1, Dim0>(tensor) * ddcHelper::get<Dim0_cov, Dim1_cov>(inv_tensor)
              + ddcHelper::get<Dim1, Dim1>(tensor) * ddcHelper::get<Dim1_cov, Dim1_cov>(inv_tensor);
    EXPECT_NEAR(id_val11, 1., TOL);
}

} // namespace

class InverseMetricTensor : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
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
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        check_inverse_tensor(metric_tensor(coords(irp)), metric_tensor.inverse(coords(irp)));
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
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        check_inverse_tensor(metric_tensor(coords(irp)), metric_tensor.inverse(coords(irp)));
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InverseMetricTensor,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(64)));
