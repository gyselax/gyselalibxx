// SPDX-License-Identifier: MIT
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
#include "toroidal_to_cylindrical.hpp"



/**
 * @brief A class for the Google tests.
 */
class InvJacobianMatrix : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};

class InvJacobianMatrix3D
    : public testing::TestWithParam<std::tuple<std::size_t, std::size_t, std::size_t>>
{
};


TEST_P(InvJacobianMatrix, InverseMatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<R, Theta, X, Y> mapping;

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_jacobian_v<CircularToCartesian<R, Theta, X, Y>, CoordRTheta>);
    InverseJacobianMatrix inv_jacobian(mapping);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        check_inverse_tensor(
                mapping.jacobian_matrix(coords(irtheta)),
                inv_jacobian(coords(irtheta)),
                1e-15);
    });
}


TEST_P(InvJacobianMatrix, InverseMatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 1.4);

    FieldMemRTheta_host<CoordRTheta> coords = get_example_coords(IdxStepR(Nr), IdxStepTheta(Nt));
    IdxRangeRTheta grid = get_idx_range(coords);

    static_assert(has_jacobian_v<CzarnyToCartesian<R, Theta, X, Y>, CoordRTheta>);
    static_assert(has_inv_jacobian_v<CzarnyToCartesian<R, Theta, X, Y>, CoordRTheta>);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        check_inverse_tensor(
                mapping.jacobian_matrix(coords(irtheta)),
                mapping.inv_jacobian_matrix(coords(irtheta)),
                1e-15);
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

    static_assert(has_jacobian_v<decltype(mapping), CoordRTheta>);
    InverseJacobianMatrix inv_jacobian(mapping);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        const CoordRTheta coord_rtheta(ddc::coordinate(irtheta));
        const double r = ddc::get<R>(coord_rtheta);
        if (fabs(r) > 1e-15) {
            check_inverse_tensor(
                    mapping.jacobian_matrix(coord_rtheta),
                    inv_jacobian(coord_rtheta),
                    1e-15);
        }
    });
}

TEST_P(InvJacobianMatrix3D, InverseMatrixToroidalDiscCzarMap)
{
    auto const [Nr, Nt, Np] = GetParam();
    const CzarnyToCartesian<R, Theta, X, Y> analytical_mapping_2d(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    CoordPhi const phi_min(0.0);
    CoordPhi const phi_max(2.0 * M_PI);
    IdxStepPhi const phi_size(Nt);

    std::vector<CoordR> r_knots = build_uniform_break_points(r_min, r_max, r_size);
    std::vector<CoordTheta> theta_knots
            = build_uniform_break_points(theta_min, theta_max, theta_size);

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(InterpPointsTheta::get_sampling<GridTheta>());
    ddc::init_discrete_space<GridPhi>(GridPhi::init<GridPhi>(phi_min, phi_max, phi_size));

    IdxRangeR interpolation_idx_range_r(InterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_theta(InterpPointsTheta::get_domain<GridTheta>());
    IdxRangePhi interpolation_idx_range_phi(IdxPhi(0), phi_size);
    IdxRangeRTheta idx_range_rtheta(interpolation_idx_range_r, interpolation_idx_range_theta);
    IdxRangeRThetaPhi idx_range(
            interpolation_idx_range_r,
            interpolation_idx_range_theta,
            interpolation_idx_range_phi);

    SplineRThetaBuilder_host builder(idx_range_rtheta);
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
                    analytical_mapping_2d,
                    builder,
                    evaluator);

    using Mapping2D = DiscreteToCartesian<
            X,
            Y,
            SplineRThetaEvaluator_host,
            R,
            Theta,
            typename Kokkos::DefaultHostExecutionSpace::memory_space>;
    DiscreteToCartesian mapping_2d = mapping_builder();

    ToroidalToCylindrical<Mapping2D, Zeta, Phi> mapping(mapping_2d);

    static_assert(has_2d_jacobian_v<decltype(mapping), CoordRThetaPhi>);
    InverseJacobianMatrix inv_jacobian(mapping);

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(idx_range, [&](IdxRThetaPhi const idx) {
        const CoordRThetaPhi coord(ddc::coordinate(idx));
        const double r = ddc::get<R>(coord);
        if (fabs(r) > 1e-15) {
            check_inverse_tensor(mapping.jacobian_matrix(coord), inv_jacobian(coord), 1e-15);
        }
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InvJacobianMatrix,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(128)));


INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InvJacobianMatrix3D,
        testing::Combine(
                testing::Values<std::size_t>(32),
                testing::Values<std::size_t>(32),
                testing::Values<std::size_t>(8)));
