// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "geometry_r_theta.hpp"

// --- Spline definitions
int constexpr BSDegreeR = 3;
int constexpr BSDegreeTheta = 3;

bool constexpr BsplineOnUniformCellsR = false;
bool constexpr BsplineOnUniformCellsTheta = false;

struct BSplinesR
    : std::conditional_t<
              BsplineOnUniformCellsR,
              ddc::UniformBSplines<R, BSDegreeR>,
              ddc::NonUniformBSplines<R, BSDegreeR>>
{
};
struct BSplinesTheta
    : std::conditional_t<
              BsplineOnUniformCellsTheta,
              ddc::UniformBSplines<Theta, BSDegreeTheta>,
              ddc::NonUniformBSplines<Theta, BSDegreeTheta>>
{
};
struct PolarBSplinesRTheta : PolarBSplines<BSplinesR, BSplinesTheta, 1>
{
};

ddc::SplineBuilderClosure constexpr SplineRClosure = ddc::SplineBuilderClosure::GREVILLE;
ddc::SplineBuilderClosure constexpr SplineThetaClosure = ddc::SplineBuilderClosure::PERIODIC;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRClosure, SplineRClosure>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaClosure, SplineThetaClosure>;

// --- Operators
using SplineRThetaBuilder_host = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::HostSpace,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        SplineRClosure, // boundary at r=0
        SplineRClosure, // boundary at rmax
        SplineThetaClosure,
        SplineThetaClosure,
        ddc::SplineSolver::LAPACK>;

using SplineRThetaEvaluatorConstBound_host = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::HostSpace,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at r=0
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>>;

using SplineRThetaEvaluatorNullBound_host = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::HostSpace,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule, // boundary at r=0
        ddc::NullExtrapolationRule, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>>;

using SplineRThetaBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        SplineRClosure, // boundary at r=0
        SplineRClosure, // boundary at rmax
        SplineThetaClosure,
        SplineThetaClosure,
        ddc::SplineSolver::LAPACK>;

using SplineRThetaEvaluatorConstBound = ddc::SplineEvaluator2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at r=0
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>>;

using SplineRThetaEvaluatorNullBound = ddc::SplineEvaluator2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule, // boundary at r=0
        ddc::NullExtrapolationRule, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>>;

using IdxRangeBSR = IdxRange<BSplinesR>;
using IdxRangeBSTheta = IdxRange<BSplinesTheta>;
using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;
using IdxRangeBSPolar = IdxRange<PolarBSplinesRTheta>;

// --- Spline representation definitions
using Spline2DMem = DFieldMem<IdxRangeBSRTheta>;
using Spline2D = DField<IdxRangeBSRTheta>;
using ConstSpline2D = DConstField<IdxRangeBSRTheta>;

/**
 * @brief Tag the polar B-splines decomposition of a function.
 *
 * Store the polar B-splines coefficients of the function.
 */
using PolarSplineMemRTheta = DFieldMem<IdxRange<PolarBSplinesRTheta>>;
using PolarSplineRTheta = DField<IdxRange<PolarBSplinesRTheta>>;

/**
 * @brief Type of the index of an element of polar B-splines.
 */
using IdxPolarBspl = Idx<PolarBSplinesRTheta>;



template <class Dim1, class Dim2>
using VectorSplineCoeffsMem2D
        = VectorFieldMem<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using VectorSplineCoeffs2D = VectorField<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using ConstVectorSplineCoeffs2D
        = VectorConstField<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2>>;
