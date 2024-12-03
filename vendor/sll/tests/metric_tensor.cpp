/// Test of the metric tensor and its inverse: (singular point avoided)
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/metric_tensor.hpp>
#include <sll/polar_bsplines.hpp>
#include <sll/view.hpp>

#include "sll/mapping/circular_to_cartesian.hpp"
#include "sll/mapping/czarny_to_cartesian.hpp"

#include "test_utils.hpp"

struct X
{
    static bool constexpr PERIODIC = false;
};
struct Y
{
    static bool constexpr PERIODIC = false;
};
struct R
{
    static bool constexpr PERIODIC = false;
};

struct Theta
{
    static bool constexpr PERIODIC = true;
};

using CoordR = ddc::Coordinate<R>;
using CoordTheta = ddc::Coordinate<Theta>;
using CoordRTheta = ddc::Coordinate<R, Theta>;

int constexpr BSDegree = 3;

struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
{
};
struct BSplinesTheta : ddc::NonUniformBSplines<Theta, BSDegree>
{
};
struct PolarBSplinesRTheta : PolarBSplines<BSplinesR, BSplinesTheta, 1>
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

using BSIdxRangeR = ddc::DiscreteDomain<BSplinesR>;
using BSIdxRangeTheta = ddc::DiscreteDomain<BSplinesTheta>;
using BSIdxRangeRTheta = ddc::DiscreteDomain<BSplinesR, BSplinesTheta>;
using BSIdxRangePolar = ddc::DiscreteDomain<PolarBSplinesRTheta>;

using IdxR = ddc::DiscreteElement<GridR>;
using IdxTheta = ddc::DiscreteElement<GridTheta>;
using IdxRTheta = ddc::DiscreteElement<GridR, GridTheta>;

using IdxStepR = ddc::DiscreteVector<GridR>;
using IdxStepTheta = ddc::DiscreteVector<GridTheta>;
using IdxStepRTheta = ddc::DiscreteVector<GridR, GridTheta>;

using IdxRangeRTheta = ddc::DiscreteDomain<GridR, GridTheta>;


template <class ElementType>
using FieldMemRTheta = ddc::Chunk<ElementType, IdxRangeRTheta>;


using Matrix_2x2 = std::array<std::array<double, 2>, 2>;


namespace {

void check_inverse(Matrix_2x2 matrix, Matrix_2x2 inv)
{
    double TOL = 1e-10;
    std::size_t N = 2;

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            double id_val = 0.0;
            for (std::size_t k(0); k < N; ++k) {
                id_val += matrix[i][k] * inv[j][k];
            }
            EXPECT_NEAR(id_val, static_cast<double>(i == j), TOL);
        }
    }
}

} // namespace

class InverseMetricTensor : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};

TEST_P(InverseMetricTensor, InverseMatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    IdxR const r_start(1); // avoid singular point.
    IdxTheta const theta_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((theta_max - theta_min) / theta_size);

    ddc::DiscreteDomain<GridR> idx_range_r(r_start, r_size);
    ddc::DiscreteDomain<GridTheta> idx_range_theta(theta_start, theta_size);
    ddc::DiscreteDomain<GridR, GridTheta> grid(idx_range_r, idx_range_theta);

    FieldMemRTheta<CoordRTheta> coords(grid);
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        coords(irp) = CoordRTheta(
                r_min + dr * ddc::select<GridR>(irp).uid(),
                theta_min + dp * ddc::select<GridR>(irp).uid());
    });

    MetricTensor<CircularToCartesian<R, Theta, X, Y>, CoordRTheta> metric_tensor(mapping);
    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        Matrix_2x2 matrix;
        Matrix_2x2 inv_matrix;

        metric_tensor(matrix, coords(irp));
        metric_tensor.inverse(inv_matrix, coords(irp));

        check_inverse(matrix, inv_matrix);
    });
}



TEST_P(InverseMetricTensor, InverseMatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    IdxR const r_start(1); // avoid singular point.
    IdxTheta const theta_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((theta_max - theta_min) / theta_size);

    ddc::DiscreteDomain<GridR> idx_range_r(r_start, r_size);
    ddc::DiscreteDomain<GridTheta> idx_range_theta(theta_start, theta_size);
    ddc::DiscreteDomain<GridR, GridTheta> grid(idx_range_r, idx_range_theta);

    FieldMemRTheta<CoordRTheta> coords(grid);
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        coords(irp) = CoordRTheta(
                r_min + dr * ddc::select<GridR>(irp).uid(),
                theta_min + dp * ddc::select<GridR>(irp).uid());
    });

    MetricTensor<CzarnyToCartesian<R, Theta, X, Y>, CoordRTheta> metric_tensor(mapping);
    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        Matrix_2x2 matrix;
        Matrix_2x2 inv_matrix;

        metric_tensor(matrix, coords(irp));
        metric_tensor.inverse(inv_matrix, coords(irp));

        check_inverse(matrix, inv_matrix);
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InverseMetricTensor,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(64)));
