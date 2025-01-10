#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_builder.hpp>
#include <sll/mapping/discrete_to_cartesian.hpp>
#include <sll/mapping/inverse_jacobian_matrix.hpp>

#include <gtest/gtest.h>



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

using SplineRThetaBuilder_host = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridR,
        GridTheta>;

using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule,
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        GridR,
        GridTheta>;

using BSIdxRangeR = ddc::DiscreteDomain<BSplinesR>;
using BSIdxRangeTheta = ddc::DiscreteDomain<BSplinesTheta>;
using BSIdxRangeRTheta = ddc::DiscreteDomain<BSplinesR, BSplinesTheta>;

using IdxRangeR = ddc::DiscreteDomain<GridR>;
using IdxRangeTheta = ddc::DiscreteDomain<GridTheta>;
using IdxRangeRTheta = ddc::DiscreteDomain<GridR, GridTheta>;

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

/**
 * @brief Check if the two matrix given as input are the same.
 *
 * The error tolerance is given at 1e-15.
 *
 * @param[in] matrix
 * 			The matrix defined with the matrix function.
 * @param[in] matrix_coeff
 * 			The matrix defined with the matrix coefficient functions.
 */
void check_matrix(Matrix_2x2 matrix, Matrix_2x2 matrix_coeff)
{
    const double TOL = 1e-15;
    constexpr std::size_t N = 2;

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            const double id_val = fabs(matrix[i][j] - matrix_coeff[i][j]);
            EXPECT_NEAR(id_val, 0., TOL);
        }
    }
}

} // namespace


/**
 * @brief A class for the Google tests.
 */
class JacobianMatrixAndJacobianCoefficients
    : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};


// Circular mapping ------------------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    IdxR const r_start(1); // avoid singular point at r = 0.
    IdxTheta const theta_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dtheta((theta_max - theta_min) / theta_size);

    ddc::DiscreteDomain<GridR> idx_range_r(r_start, r_size);
    ddc::DiscreteDomain<GridTheta> idx_range_theta(theta_start, theta_size);
    ddc::DiscreteDomain<GridR, GridTheta> grid(idx_range_r, idx_range_theta);

    FieldMemRTheta<CoordRTheta> coords(grid);
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        coords(irp) = CoordRTheta(
                r_min + dr * ddc::select<GridR>(irp).uid(),
                theta_min + dtheta * ddc::select<GridR>(irp).uid());
    });

    static_assert(has_2d_jacobian_v<decltype(mapping), CoordRTheta>);
    static_assert(has_2d_inv_jacobian_v<decltype(mapping), CoordRTheta>);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 Jacobian_matrix_coeff;

        mapping.jacobian_matrix(coords(irp), Jacobian_matrix);
        Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coords(irp));
        Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coords(irp));
        Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coords(irp));
        Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coords(irp));

        check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);
    });

    // --- for the inverse Jacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        Matrix_2x2 inv_Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix_coeff;

        mapping.inv_jacobian_matrix(coords(irp), inv_Jacobian_matrix);
        inv_Jacobian_matrix_coeff[0][0] = mapping.inv_jacobian_11(coords(irp));
        inv_Jacobian_matrix_coeff[0][1] = mapping.inv_jacobian_12(coords(irp));
        inv_Jacobian_matrix_coeff[1][0] = mapping.inv_jacobian_21(coords(irp));
        inv_Jacobian_matrix_coeff[1][1] = mapping.inv_jacobian_22(coords(irp));

        check_matrix(inv_Jacobian_matrix, inv_Jacobian_matrix_coeff);
    });
}



// Czarny mapping --------------------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    IdxR const r_start(1); // avoid singular point at r = 0.
    IdxTheta const theta_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dtheta((theta_max - theta_min) / theta_size);

    ddc::DiscreteDomain<GridR> idx_range_r(r_start, r_size);
    ddc::DiscreteDomain<GridTheta> idx_range_theta(theta_start, theta_size);
    ddc::DiscreteDomain<GridR, GridTheta> grid(idx_range_r, idx_range_theta);

    FieldMemRTheta<CoordRTheta> coords(grid);
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        coords(irp) = CoordRTheta(
                r_min + dr * ddc::select<GridR>(irp).uid(),
                theta_min + dtheta * ddc::select<GridR>(irp).uid());
    });

    static_assert(has_2d_jacobian_v<decltype(mapping), CoordRTheta>);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 Jacobian_matrix_coeff;

        mapping.jacobian_matrix(coords(irp), Jacobian_matrix);
        Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coords(irp));
        Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coords(irp));
        Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coords(irp));
        Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coords(irp));

        check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);
    });

    InverseJacobianMatrix inv_jacobian(mapping);
    // --- for the inverseJacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        Matrix_2x2 inv_Jacobian_matrix = inv_jacobian(coords(irp));
        Matrix_2x2 inv_Jacobian_matrix_coeff;

        inv_Jacobian_matrix_coeff[0][0] = inv_jacobian.inv_jacobian_11(coords(irp));
        inv_Jacobian_matrix_coeff[0][1] = inv_jacobian.inv_jacobian_12(coords(irp));
        inv_Jacobian_matrix_coeff[1][0] = inv_jacobian.inv_jacobian_21(coords(irp));
        inv_Jacobian_matrix_coeff[1][1] = inv_jacobian.inv_jacobian_22(coords(irp));

        check_matrix(inv_Jacobian_matrix, inv_Jacobian_matrix_coeff);
    });
}



// Discrete Czarny mapping -----------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixDiscCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> analytical_mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dtheta((theta_max - theta_min) / theta_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordTheta> theta_knots(theta_size + 1);

    for (int i(0); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    r_knots[r_size] = CoordR(r_max);
    for (int i(0); i < theta_size + 1; ++i) {
        theta_knots[i] = CoordTheta(theta_min + i * dtheta);
    }

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(InterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR interpolation_idx_range_R(InterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_Theta(InterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_R, interpolation_idx_range_Theta);

    SplineRThetaBuilder_host builder(grid);
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluator evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);
    DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator>
            mapping_builder(
                    Kokkos::DefaultHostExecutionSpace(),
                    analytical_mapping,
                    builder,
                    evaluator);
    DiscreteToCartesian mapping = mapping_builder();

    static_assert(has_2d_jacobian_v<decltype(mapping), CoordRTheta>);
    InverseJacobianMatrix inv_jacobian(mapping);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        const CoordRTheta coord_rp(ddc::coordinate(irp));
        const double r = ddc::get<R>(coord_rp);
        if (fabs(r) > 1e-15) {
            // --- for the Jacobian matrix:
            Matrix_2x2 Jacobian_matrix;
            Matrix_2x2 Jacobian_matrix_coeff;

            mapping.jacobian_matrix(coord_rp, Jacobian_matrix);
            Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coord_rp);
            Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coord_rp);
            Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coord_rp);
            Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coord_rp);

            check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);


            // --- for the inverse Jacobian matrix:
            Matrix_2x2 inv_Jacobian_matrix = inv_jacobian(coord_rp);
            Matrix_2x2 inv_Jacobian_matrix_coeff;

            inv_Jacobian_matrix_coeff[0][0] = inv_jacobian.inv_jacobian_11(coord_rp);
            inv_Jacobian_matrix_coeff[0][1] = inv_jacobian.inv_jacobian_12(coord_rp);
            inv_Jacobian_matrix_coeff[1][0] = inv_jacobian.inv_jacobian_21(coord_rp);
            inv_Jacobian_matrix_coeff[1][1] = inv_jacobian.inv_jacobian_22(coord_rp);

            check_matrix(inv_Jacobian_matrix, inv_Jacobian_matrix_coeff);
        }
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        JacobianMatrixAndJacobianCoefficients,
        testing::Combine(testing::Values<std::size_t>(40), testing::Values<std::size_t>(80)));
