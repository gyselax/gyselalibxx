/// Test of the metric tensor and its inverse: (singular point avoided)
#include <sll/bsplines_non_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/polar_bsplines.hpp>

#include "sll/mapping/circular_to_cartesian.hpp"
#include "sll/mapping/czarny_to_cartesian.hpp"

#include "test_utils.hpp"

struct DimX
{
    static bool constexpr PERIODIC = false;
};
struct DimY
{
    static bool constexpr PERIODIC = false;
};
struct DimR
{
    static bool constexpr PERIODIC = false;
};

struct DimP
{
    static bool constexpr PERIODIC = true;
};

using CoordR = Coordinate<DimR>;
using CoordP = Coordinate<DimP>;
using CoordRP = Coordinate<DimR, DimP>;

int constexpr BSDegree = 3;

using BSplinesR = NonUniformBSplines<DimR, BSDegree>;
using BSplinesP = NonUniformBSplines<DimP, BSDegree>;
using PolarBSplinesRP = PolarBSplines<BSplinesR, BSplinesP, 1>;

using InterpPointsR
        = GrevilleInterpolationPoints<BSplinesR, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using InterpPointsP
        = GrevilleInterpolationPoints<BSplinesP, BoundCond::PERIODIC, BoundCond::PERIODIC>;

using IDimR = typename InterpPointsR::interpolation_mesh_type;
using IDimP = typename InterpPointsP::interpolation_mesh_type;

using BSDomainR = DiscreteDomain<BSplinesR>;
using BSDomainP = DiscreteDomain<BSplinesP>;
using BSDomainRP = DiscreteDomain<BSplinesR, BSplinesP>;
using BSDomainPolar = DiscreteDomain<PolarBSplinesRP>;

using IndexR = DiscreteElement<IDimR>;
using IndexP = DiscreteElement<IDimP>;
using IndexRP = DiscreteElement<IDimR, IDimP>;

using IVectR = DiscreteVector<IDimR>;
using IVectP = DiscreteVector<IDimP>;
using IVectRP = DiscreteVector<IDimR, IDimP>;

using IDomainRP = DiscreteDomain<IDimR, IDimP>;


template <class ElementType>
using FieldRP = Chunk<ElementType, IDomainRP>;


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
    const CircularToCartesian<DimX, DimY, DimR, DimP> mapping;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    IndexR const r_start(1); // avoid singular point.
    IndexP const p_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    ddc::DiscreteDomain<IDimR> domain_r(r_start, r_size);
    ddc::DiscreteDomain<IDimP> domain_p(p_start, p_size);
    ddc::DiscreteDomain<IDimR, IDimP> grid(domain_r, domain_p);

    FieldRP<CoordRP> coords(grid);
    for_each(grid, [&](IndexRP const irp) {
        coords(irp) = CoordRP(
                r_min + dr * select<IDimR>(irp).uid(),
                p_min + dp * select<IDimR>(irp).uid());
    });

    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    for_each(grid, [&](IndexRP const irp) {
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
    const CzarnyToCartesian<DimX, DimY, DimR, DimP> mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    IndexR const r_start(1); // avoid singular point.
    IndexP const p_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    ddc::DiscreteDomain<IDimR> domain_r(r_start, r_size);
    ddc::DiscreteDomain<IDimP> domain_p(p_start, p_size);
    ddc::DiscreteDomain<IDimR, IDimP> grid(domain_r, domain_p);

    FieldRP<CoordRP> coords(grid);
    for_each(grid, [&](IndexRP const irp) {
        coords(irp) = CoordRP(
                r_min + dr * select<IDimR>(irp).uid(),
                p_min + dp * select<IDimR>(irp).uid());
    });

    // Test for each coordinates if the inverse_metric_tensor is the inverse of the metric_tensor
    for_each(grid, [&](IndexRP const irp) {
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
