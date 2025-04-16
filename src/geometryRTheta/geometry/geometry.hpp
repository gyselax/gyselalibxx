// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "polar_bsplines.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"


/*
 * @file geometry.hpp
 *
 * Definition of
 *   - @f$ r@f$, @f$ \theta@f$, @f$(r, \theta)@f$ dimensions.
 *   - @f$x@f$, @f$y@f$, @f$(x, y)@f$ dimensions.
 */


// POLAR SPACE AND VELOCITY ----------------------------------------------------------------------
// --- Continuous dimensions
struct R_cov;
struct Theta_cov;
/**
 * @brief Define non periodic real contravariant R dimension.
 */
struct R
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = false;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = R_cov;
};
/**
 * @brief Define periodic real contravariant Theta dimension.
 */
struct Theta
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = false;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = Theta_cov;
};

/**
 * @brief Define non periodic real covariant R dimension.
 */
struct R_cov
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = false;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = R;
};
/**
 * @brief Define periodic real covariant Theta dimension.
 */
struct Theta_cov
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = false;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = Theta;
};

/**
 * @brief Define non periodic real R velocity dimension.
 */
struct Vr
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define periodic real Theta velocity dimension.
 */
struct Vtheta
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;
};


using CoordR = Coord<R>;
using CoordTheta = Coord<Theta>;
using CoordRTheta = Coord<R, Theta>;

using CoordVr = Coord<Vr>;
using CoordVtheta = Coord<Vtheta>;

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

ddc::BoundCond constexpr SplineRBoundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineThetaBoundary = ddc::BoundCond::PERIODIC;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaBoundary, SplineThetaBoundary>;

// --- Discrete dimensions
struct GridR : SplineInterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : SplineInterpPointsTheta::interpolation_discrete_dimension_type
{
};

// --- Operators
using SplineRThetaBuilder_host = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::HostSpace,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        SplineRBoundary, // boundary at r=0
        SplineRBoundary, // boundary at rmax
        SplineThetaBoundary,
        SplineThetaBoundary,
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
        SplineRBoundary, // boundary at r=0
        SplineRBoundary, // boundary at rmax
        SplineThetaBoundary,
        SplineThetaBoundary,
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
// --- Index definitions
using IdxR = Idx<GridR>;
using IdxTheta = Idx<GridTheta>;
using IdxRTheta = Idx<GridR, GridTheta>;

// --- Index Step definitions
using IdxStepR = IdxStep<GridR>;
using IdxStepTheta = IdxStep<GridTheta>;
using IdxStepRTheta = IdxStep<GridR, GridTheta>;

// --- Index Range definitions
using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;

using IdxRangeBSR = IdxRange<BSplinesR>;
using IdxRangeBSTheta = IdxRange<BSplinesTheta>;
using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;
using IdxRangeBSPolar = IdxRange<PolarBSplinesRTheta>;


// --- FieldMem definitions
template <class ElementType>
using FieldMemR = FieldMem<ElementType, IdxRangeR>;

template <class ElementType>
using FieldMemTheta = FieldMem<ElementType, IdxRangeTheta>;

template <class ElementType>
using FieldMemRTheta = FieldMem<ElementType, IdxRangeRTheta>;

using DFieldMemR = FieldMemR<double>;
using DFieldMemTheta = FieldMemTheta<double>;
using DFieldMemRTheta = FieldMemRTheta<double>;

// --- Field definitions
template <class ElementType>
using FieldR = Field<ElementType, IdxRangeR>;

template <class ElementType>
using FieldTheta = Field<ElementType, IdxRangeTheta>;

template <class ElementType>
using FieldRTheta = Field<ElementType, IdxRangeRTheta>;

using DFieldR = FieldR<double>;
using DFieldTheta = FieldTheta<double>;
using DFieldRTheta = FieldRTheta<double>;

// --- Const Field definitions
template <class ElementType>
using ConstFieldR = ConstField<ElementType, IdxRangeR>;

template <class ElementType>
using ConstFieldTheta = ConstField<ElementType, IdxRangeTheta>;

template <class ElementType>
using ConstFieldRTheta = ConstField<ElementType, IdxRangeRTheta>;

using DConstFieldR = ConstFieldR<double>;
using DConstFieldTheta = ConstFieldTheta<double>;
using DConstFieldRTheta = ConstFieldRTheta<double>;

// --- Spline representation definitions
using Spline2DMem = DFieldMem<IdxRangeBSRTheta>;
using Spline2D = DField<IdxRangeBSRTheta>;
using ConstSpline2D = DConstField<IdxRangeBSRTheta>;

/**
 * @brief Tag the polar B-splines decomposition of a function.
 *
 * Store the polar B-splines coefficients of the function.
 */
using PolarSplineMemRTheta = PolarSplineMem<PolarBSplinesRTheta>;

/**
 * @brief Type of the index of an element of polar B-splines.
 */
using IdxPolarBspl = Idx<PolarBSplinesRTheta>;


// --- VectorFieldMem definitions
template <class Dim1, class Dim2>
using DVectorFieldMemRTheta = VectorFieldMem<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using DVectorFieldRTheta = VectorField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using DConstVectorFieldRTheta
        = VectorConstField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>>;



template <class Dim1, class Dim2>
using VectorSplineCoeffsMem2D
        = VectorFieldMem<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using VectorSplineCoeffs2D = VectorField<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using ConstVectorSplineCoeffs2D
        = VectorConstField<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2>>;



// CARTESIAN SPACE AND VELOCITY ------------------------------------------------------------------
// --- Continuous dimensions
/**
 * @brief Define non periodic real X dimension.
 */
struct X
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = X;
};
/**
 * @brief Define non periodic real Y dimension.
 */
struct Y
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = Y;
};

/**
 * @brief Define non periodic real X velocity dimension.
 */
struct Vx
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define non periodic real Y velocity dimension.
 */
struct Vy
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};


using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

using CoordVx = Coord<Vx>;
using CoordVy = Coord<Vy>;
