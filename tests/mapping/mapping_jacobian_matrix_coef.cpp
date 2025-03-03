#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "mapping_testing_tools.hpp"
#include "mesh_builder.hpp"



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

using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;

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

/**
 * @brief Check if the two matrix given as input are the same.
 *
 * The error tolerance is given at 1e-14.
 *
 * @param[in] matrix
 * 			The matrix defined with the matrix function.
 * @param[in] matrix_coeff
 * 			The matrix defined with the matrix coefficient functions.
 */
void check_matrix(Matrix_2x2 matrix, Matrix_2x2 matrix_coeff)
{
    const double TOL = 1e-13;
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

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_2d_jacobian_v<decltype(mapping), CoordRTheta>);
    static_assert(has_2d_inv_jacobian_v<decltype(mapping), CoordRTheta>);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 Jacobian_matrix_coeff;

        mapping.jacobian_matrix(coords(irtheta), Jacobian_matrix);
        Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coords(irtheta));
        Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coords(irtheta));
        Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coords(irtheta));
        Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coords(irtheta));

        check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);
    });

    // --- for the inverse Jacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        Matrix_2x2 inv_Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix_coeff;

        mapping.inv_jacobian_matrix(coords(irtheta), inv_Jacobian_matrix);
        inv_Jacobian_matrix_coeff[0][0] = mapping.inv_jacobian_11(coords(irtheta));
        inv_Jacobian_matrix_coeff[0][1] = mapping.inv_jacobian_12(coords(irtheta));
        inv_Jacobian_matrix_coeff[1][0] = mapping.inv_jacobian_21(coords(irtheta));
        inv_Jacobian_matrix_coeff[1][1] = mapping.inv_jacobian_22(coords(irtheta));

        check_matrix(inv_Jacobian_matrix, inv_Jacobian_matrix_coeff);
    });
}



// Czarny mapping --------------------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 1.4);

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_2d_jacobian_v<decltype(mapping), CoordRTheta>);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 Jacobian_matrix_coeff;

        mapping.jacobian_matrix(coords(irtheta), Jacobian_matrix);
        Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coords(irtheta));
        Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coords(irtheta));
        Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coords(irtheta));
        Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coords(irtheta));

        check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);
    });

    InverseJacobianMatrix inv_jacobian(mapping);
    // --- for the inverseJacobian matrix:
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        Matrix_2x2 inv_Jacobian_matrix = inv_jacobian(coords(irtheta));
        Matrix_2x2 inv_Jacobian_matrix_coeff;

        inv_Jacobian_matrix_coeff[0][0] = inv_jacobian.inv_jacobian_11(coords(irtheta));
        inv_Jacobian_matrix_coeff[0][1] = inv_jacobian.inv_jacobian_12(coords(irtheta));
        inv_Jacobian_matrix_coeff[1][0] = inv_jacobian.inv_jacobian_21(coords(irtheta));
        inv_Jacobian_matrix_coeff[1][1] = inv_jacobian.inv_jacobian_22(coords(irtheta));

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

    std::vector<CoordR> r_knots = build_uniform_break_points(r_min, r_max, r_size);
    std::vector<CoordTheta> theta_knots
            = build_uniform_break_points(theta_min, theta_max, theta_size);

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(InterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR interpolation_idx_range_R(InterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_theta(InterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_R, interpolation_idx_range_theta);

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
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        const CoordRTheta coord_rtheta(ddc::coordinate(irtheta));
        const double r = ddc::get<R>(coord_rtheta);
        if (fabs(r) > 1e-15) {
            // --- for the Jacobian matrix:
            Matrix_2x2 Jacobian_matrix;
            Matrix_2x2 Jacobian_matrix_coeff;

            mapping.jacobian_matrix(coord_rtheta, Jacobian_matrix);
            Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coord_rtheta);
            Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coord_rtheta);
            Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coord_rtheta);
            Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coord_rtheta);

            check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);


            // --- for the inverse Jacobian matrix:
            Matrix_2x2 inv_Jacobian_matrix = inv_jacobian(coord_rtheta);
            Matrix_2x2 inv_Jacobian_matrix_coeff;

            inv_Jacobian_matrix_coeff[0][0] = inv_jacobian.inv_jacobian_11(coord_rtheta);
            inv_Jacobian_matrix_coeff[0][1] = inv_jacobian.inv_jacobian_12(coord_rtheta);
            inv_Jacobian_matrix_coeff[1][0] = inv_jacobian.inv_jacobian_21(coord_rtheta);
            inv_Jacobian_matrix_coeff[1][1] = inv_jacobian.inv_jacobian_22(coord_rtheta);

            check_matrix(inv_Jacobian_matrix, inv_Jacobian_matrix_coeff);
        }
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        JacobianMatrixAndJacobianCoefficients,
        testing::Combine(testing::Values<std::size_t>(40), testing::Values<std::size_t>(80)));
