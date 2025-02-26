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
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = X;
};
struct Y
{
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
 * @brief Check if the product of the matrix and inv_matrix gives the identity matrix.
 *
 * The error tolerance is given at 1e-15.
 *
 * @param[in] matrix
 * 			The Jacobian matrix of the mapping.
 * @param[in] inv_matrix
 * 			The inverse Jacobian matrix of the mapping.
 */
template <class StartDims, class EndDims>
void check_inverse_tensor(
        DTensor<StartDims, EndDims> const& tensor,
        DTensor<vector_index_set_dual_t<EndDims>, vector_index_set_dual_t<StartDims>> const&
                inv_tensor)
{
    double TOL = 1e-10;

    using StartDim0 = ddc::type_seq_element_t<0, StartDims>;
    using StartDim1 = ddc::type_seq_element_t<1, StartDims>;
    using StartDim0_cov = typename StartDim0::Dual;
    using StartDim1_cov = typename StartDim1::Dual;

    using EndDim0_cov = ddc::type_seq_element_t<0, EndDims>;
    using EndDim1_cov = ddc::type_seq_element_t<1, EndDims>;
    using EndDim0 = typename EndDim0_cov::Dual;
    using EndDim1 = typename EndDim1_cov::Dual;

    double const id_val00 = ddcHelper::get<StartDim0, EndDim0_cov>(tensor)
                                    * ddcHelper::get<EndDim0, StartDim0_cov>(inv_tensor)
                            + ddcHelper::get<StartDim0, EndDim1_cov>(tensor)
                                      * ddcHelper::get<EndDim1, StartDim0_cov>(inv_tensor);
    EXPECT_NEAR(id_val00, 1., TOL);

    double const id_val01 = ddcHelper::get<StartDim0, EndDim0_cov>(tensor)
                                    * ddcHelper::get<EndDim0, StartDim1_cov>(inv_tensor)
                            + ddcHelper::get<StartDim0, EndDim1_cov>(tensor)
                                      * ddcHelper::get<EndDim1, StartDim1_cov>(inv_tensor);
    EXPECT_NEAR(id_val01, 0., TOL);

    double const id_val10 = ddcHelper::get<StartDim1, EndDim0_cov>(tensor)
                                    * ddcHelper::get<EndDim0, StartDim0_cov>(inv_tensor)
                            + ddcHelper::get<StartDim1, EndDim1_cov>(tensor)
                                      * ddcHelper::get<EndDim1, StartDim0_cov>(inv_tensor);
    EXPECT_NEAR(id_val10, 0., TOL);

    double const id_val11 = ddcHelper::get<StartDim1, EndDim0_cov>(tensor)
                                    * ddcHelper::get<EndDim0, StartDim1_cov>(inv_tensor)
                            + ddcHelper::get<StartDim1, EndDim1_cov>(tensor)
                                      * ddcHelper::get<EndDim1, StartDim1_cov>(inv_tensor);
    EXPECT_NEAR(id_val11, 1., TOL);
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
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        check_inverse_tensor(mapping.jacobian_matrix(coords(irp)), inv_jacobian(coords(irp)));
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
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        check_inverse_tensor(
                mapping.jacobian_matrix(coords(irp)),
                mapping.inv_jacobian_matrix(coords(irp)));
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

    // Test for each coordinates if the inv_Jacobian_matrix is the inverse of the Jacobian_matrix
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        const CoordRTheta coord_rp(ddc::coordinate(irp));
        const double r = ddc::get<R>(coord_rp);
        if (fabs(r) > 1e-15) {
            check_inverse_tensor(mapping.jacobian_matrix(coord_rp), inv_jacobian(coord_rp));
        }
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        InvJacobianMatrix,
        testing::Combine(testing::Values<std::size_t>(64), testing::Values<std::size_t>(128)));
