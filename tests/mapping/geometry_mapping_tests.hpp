// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "indexed_tensor.hpp"
#include "tensor.hpp"

struct X
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    /// The corresponding type in the dual space.
    using Dual = X;
};
struct Y
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    /// The corresponding type in the dual space.
    using Dual = Y;
};
struct R_cov;
struct Theta_cov;
struct R
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    /// The corresponding type in the dual space.
    using Dual = R_cov;
};

struct Theta
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    /// The corresponding type in the dual space.
    using Dual = Theta_cov;
};
struct R_cov
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    /// The corresponding type in the dual space.
    using Dual = R;
};

struct Theta_cov
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    /// The corresponding type in the dual space.
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

template <class ExecSpace>
using SplineRThetaBuilder = ddc::SplineBuilder2D<
        ExecSpace,
        typename ExecSpace::memory_space,
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
using SplineRThetaBuilder_host = SplineRThetaBuilder<Kokkos::DefaultHostExecutionSpace>;

template <class ExecSpace>
using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        ExecSpace,
        typename ExecSpace::memory_space,
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
using SplineRThetaEvaluator_host = SplineRThetaEvaluator<Kokkos::DefaultHostExecutionSpace>;

using IdxRangeBSR = IdxRange<BSplinesR>;
using IdxRangeBSTheta = IdxRange<BSplinesTheta>;
using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;

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
 * @param[in] matrix
 * 			The Jacobian matrix of the mapping.
 * @param[in] inv_matrix
 * 			The inverse Jacobian matrix of the mapping.
 * @param[in] TOL
 *          The error tolerance.
 */
template <class StartDims, class EndDims>
void check_inverse_tensor(
        DTensor<StartDims, EndDims> const& tensor,
        DTensor<vector_index_set_dual_t<EndDims>, vector_index_set_dual_t<StartDims>> const&
                inv_tensor,
        double TOL)
{
    using StartDim0 = ddc::type_seq_element_t<0, StartDims>;
    using StartDim1 = ddc::type_seq_element_t<1, StartDims>;
    using StartDim0_cov = typename StartDim0::Dual;
    using StartDim1_cov = typename StartDim1::Dual;

    DTensor<StartDims, vector_index_set_dual_t<StartDims>> identity
            = tensor_mul(index<'i', 'j'>(tensor), index<'j', 'k'>(inv_tensor));

    EXPECT_NEAR((ddcHelper::get<StartDim0, StartDim0_cov>(identity)), 1., TOL);
    EXPECT_NEAR((ddcHelper::get<StartDim0, StartDim1_cov>(identity)), 0., TOL);
    EXPECT_NEAR((ddcHelper::get<StartDim1, StartDim0_cov>(identity)), 0., TOL);
    EXPECT_NEAR((ddcHelper::get<StartDim1, StartDim1_cov>(identity)), 1., TOL);
}
