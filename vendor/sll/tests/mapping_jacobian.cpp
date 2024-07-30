#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include "sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/circular_to_cartesian.hpp"
#include "sll/mapping/curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/czarny_to_cartesian.hpp"
#include "sll/mapping/discrete_mapping_to_cartesian.hpp"

#include "test_utils.hpp"



namespace {
struct X
{
};
struct Y
{
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

using SplineRPBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesP,
        GridR,
        GridP,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridR,
        GridP>;

using SplineRPEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesP,
        GridR,
        GridP,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule,
        ddc::PeriodicExtrapolationRule<P>,
        ddc::PeriodicExtrapolationRule<P>,
        GridR,
        GridP>;

using BSIdxRangeR = ddc::DiscreteDomain<BSplinesR>;
using BSIdxRangeP = ddc::DiscreteDomain<BSplinesP>;
using BSIdxRangeRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;

using IdxRangeR = ddc::DiscreteDomain<GridR>;
using IdxRangeP = ddc::DiscreteDomain<GridP>;
using IdxRangeRP = ddc::DiscreteDomain<GridR, GridP>;

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

/**
 * @brief Check if the product of the matrix and inv_matrix gives the identity matrix.
 *
 * The error tolerance is given at 1e-15.
 *
 * @param[in] matrix
 * 			The Jacobian matrix of the mapping.
 * @param[in] inv_matrix
 * 			The inverse Jacobian matrix of the mapping.
 */
void check_inverse(Matrix_2x2 matrix, Matrix_2x2 inv_matrix)
{
    const double TOL = 1e-15;
    constexpr std::size_t N = 2;

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            double id_val = 0.0;
            for (std::size_t k(0); k < N; ++k) {
                id_val += matrix[i][k] * inv_matrix[k][j];
            }
            EXPECT_NEAR(id_val, static_cast<double>(i == j), TOL);
        }
    }
}

} // namespace


/**
 * @brief A class for the Google tests.
 */
class InverseJacobianMatrix : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};



TEST_P(InverseJacobianMatrix, InverseMatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<X, Y, R, P> mapping;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IdxStepP const p_size(Nt);

    IdxR const r_start(1); // avoid singular point at r = 0.
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

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRP const irp) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix;

        mapping.jacobian_matrix(coords(irp), Jacobian_matrix);
        mapping.inv_jacobian_matrix(coords(irp), inv_Jacobian_matrix);

        check_inverse(Jacobian_matrix, inv_Jacobian_matrix);
    });
}


TEST_P(InverseJacobianMatrix, InverseMatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<X, Y, R, P> mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IdxStepP const p_size(Nt);

    IdxR const r_start(1); // avoid singular point at r = 0.
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
                p_min + dp * ddc::select<GridP>(irp).uid());
    });

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRP const irp) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix;

        mapping.jacobian_matrix(coords(irp), Jacobian_matrix);
        mapping.inv_jacobian_matrix(coords(irp), inv_Jacobian_matrix);

        check_inverse(Jacobian_matrix, inv_Jacobian_matrix);
    });
}


TEST_P(InverseJacobianMatrix, InverseMatrixDiscCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<X, Y, R, P> analytical_mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IdxStepP const p_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordP> p_knots(p_size + 1);

    for (int i(0); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    r_knots[r_size] = CoordR(r_max);
    for (int i(0); i < p_size + 1; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridP>(InterpPointsP::get_sampling<GridP>());

    IdxRangeR interpolation_idx_range_R(InterpPointsR::get_domain<GridR>());
    IdxRangeP interpolation_idx_range_P(InterpPointsP::get_domain<GridP>());
    IdxRangeRP grid(interpolation_idx_range_R, interpolation_idx_range_P);

    SplineRPBuilder builder(grid);
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<P> p_extrapolation_rule;
    SplineRPEvaluator evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);
    DiscreteToCartesian mapping = DiscreteToCartesian<X, Y, SplineRPBuilder, SplineRPEvaluator>::
            analytical_to_discrete(analytical_mapping, builder, evaluator);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRP const irp) {
        const CoordRP coord_rp(ddc::coordinate(irp));
        const double r = ddc::get<R>(coord_rp);
        if (fabs(r) > 1e-15) {
            Matrix_2x2 Jacobian_matrix;
            Matrix_2x2 inv_Jacobian_matrix;

            mapping.jacobian_matrix(coord_rp, Jacobian_matrix);
            mapping.inv_jacobian_matrix(coord_rp, inv_Jacobian_matrix);

            check_inverse(Jacobian_matrix, inv_Jacobian_matrix);
        }
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InverseJacobianMatrix,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(128)));
