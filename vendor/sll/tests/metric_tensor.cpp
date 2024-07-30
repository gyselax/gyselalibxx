/// Test of the metric tensor and its inverse: (singular point avoided)
#include <sll/polar_bsplines.hpp>

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

struct P
{
    static bool constexpr PERIODIC = true;
};

using CoordR = ddc::Coordinate<R>;
using CoordP = ddc::Coordinate<P>;
using CoordRP = ddc::Coordinate<R, P>;

int constexpr BSDegree = 3;

struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
{
};
struct BSplinesP : ddc::NonUniformBSplines<P, BSDegree>
{
};
struct PolarBSplinesRP : PolarBSplines<BSplinesR, BSplinesP, 1>
{
};

using InterpPointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using InterpPointsP = ddc::
        GrevilleInterpolationPoints<BSplinesP, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

struct GridR : InterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridP : InterpPointsP::interpolation_discrete_dimension_type
{
};

using BSIdxRangeR = ddc::DiscreteDomain<BSplinesR>;
using BSIdxRangeP = ddc::DiscreteDomain<BSplinesP>;
using BSIdxRangeRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;
using BSIdxRangePolar = ddc::DiscreteDomain<PolarBSplinesRP>;

using IdxR = ddc::DiscreteElement<GridR>;
using IdxP = ddc::DiscreteElement<GridP>;
using IdxRP = ddc::DiscreteElement<GridR, GridP>;

using IdxStepR = ddc::DiscreteVector<GridR>;
using IdxStepP = ddc::DiscreteVector<GridP>;
using IdxStepRP = ddc::DiscreteVector<GridR, GridP>;

using IdxRangeRP = ddc::DiscreteDomain<GridR, GridP>;


template <class ElementType>
using FieldMemRP = ddc::Chunk<ElementType, IdxRangeRP>;


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
    const CircularToCartesian<X, Y, R, P> mapping;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IdxStepP const p_size(Nt);

    IdxR const r_start(1); // avoid singular point.
    IdxP const p_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    ddc::DiscreteDomain<GridR> idx_range_r(r_start, r_size);
    ddc::DiscreteDomain<GridP> idx_range_p(p_start, p_size);
    ddc::DiscreteDomain<GridR, GridP> grid(idx_range_r, idx_range_p);

    FieldMemRP<CoordRP> coords(grid);
    ddc::for_each(grid, [&](IdxRP const irp) {
        coords(irp) = CoordRP(
                r_min + dr * ddc::select<GridR>(irp).uid(),
                p_min + dp * ddc::select<GridR>(irp).uid());
    });

    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRP const irp) {
        Matrix_2x2 matrix;
        Matrix_2x2 inv_matrix;

        mapping.metric_tensor(coords(irp), matrix);
        mapping.inverse_metric_tensor(coords(irp), inv_matrix);

        check_inverse(matrix, inv_matrix);
    });
}



TEST_P(InverseMetricTensor, InverseMatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<X, Y, R, P> mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IdxStepP const p_size(Nt);

    IdxR const r_start(1); // avoid singular point.
    IdxP const p_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    ddc::DiscreteDomain<GridR> idx_range_r(r_start, r_size);
    ddc::DiscreteDomain<GridP> idx_range_p(p_start, p_size);
    ddc::DiscreteDomain<GridR, GridP> grid(idx_range_r, idx_range_p);

    FieldMemRP<CoordRP> coords(grid);
    ddc::for_each(grid, [&](IdxRP const irp) {
        coords(irp) = CoordRP(
                r_min + dr * ddc::select<GridR>(irp).uid(),
                p_min + dp * ddc::select<GridR>(irp).uid());
    });

    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    ddc::for_each(grid, [&](IdxRP const irp) {
        Matrix_2x2 matrix;
        Matrix_2x2 inv_matrix;

        mapping.metric_tensor(coords(irp), matrix);
        mapping.inverse_metric_tensor(coords(irp), inv_matrix);

        check_inverse(matrix, inv_matrix);
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InverseMetricTensor,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(64)));
