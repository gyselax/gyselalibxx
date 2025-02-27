// SPDX-License-Identifier: MIT
#pragma once

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
using CoordXY = Coord<X, Y>;

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
