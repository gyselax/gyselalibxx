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
#include "geometry_mapping_tests.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "mapping_testing_tools.hpp"
#include "mesh_builder.hpp"


namespace {

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
class InvJacobianMatrix : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};



TEST_P(InvJacobianMatrix, InverseMatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_2d_jacobian_v<CircularToCartesian<R, Theta, X, Y>, CoordRTheta>);
    InverseJacobianMatrix inv_jacobian(mapping);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix = inv_jacobian(coords(irtheta));

        mapping.jacobian_matrix(coords(irtheta), Jacobian_matrix);

        check_inverse(Jacobian_matrix, inv_Jacobian_matrix);
    });
}


TEST_P(InvJacobianMatrix, InverseMatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 1.4);

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_2d_jacobian_v<CzarnyToCartesian<R, Theta, X, Y>, CoordRTheta>);
    static_assert(has_2d_inv_jacobian_v<CzarnyToCartesian<R, Theta, X, Y>, CoordRTheta>);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix;

        mapping.jacobian_matrix(coords(irtheta), Jacobian_matrix);
        mapping.inv_jacobian_matrix(coords(irtheta), inv_Jacobian_matrix);

        check_inverse(Jacobian_matrix, inv_Jacobian_matrix);
    });
}


TEST_P(InvJacobianMatrix, InverseMatrixDiscCzarMap)
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

    IdxRangeR interpolation_idx_range_r(InterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_theta(InterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_r, interpolation_idx_range_theta);

    SplineRThetaBuilder_host builder(grid);
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluator_host evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);
    DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator_host>
            mapping_builder(
                    Kokkos::DefaultHostExecutionSpace(),
                    analytical_mapping,
                    builder,
                    evaluator);
    DiscreteToCartesian mapping = mapping_builder();

    static_assert(has_2d_jacobian_v<decltype(mapping), CoordRTheta>);
    InverseJacobianMatrix inv_jacobian(mapping);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        const CoordRTheta coord_rtheta(ddc::coordinate(irtheta));
        const double r = ddc::get<R>(coord_rtheta);
        if (fabs(r) > 1e-15) {
            Matrix_2x2 Jacobian_matrix;
            Matrix_2x2 inv_Jacobian_matrix = inv_jacobian(coord_rtheta);

            mapping.jacobian_matrix(coord_rtheta, Jacobian_matrix);

            check_inverse(Jacobian_matrix, inv_Jacobian_matrix);
        }
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InvJacobianMatrix,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(128)));
