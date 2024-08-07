#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/polar_bsplines.hpp>

#include <directional_tag.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

#include "ddc_aliases.hpp"


/*
 * @file geometry.hpp
 *
 * Definition of
 *   - @f$ r@f$, @f$ \theta@f$, @f$(r, \theta)@f$ dimensions.
 *   - @f$x@f$, @f$y@f$, @f$(x, y)@f$ dimensions.
 */


// POLAR SPACE AND VELOCITY ----------------------------------------------------------------------
// --- Continuous dimensions
/**
 * @brief Define non periodic real R dimension.
 */
struct R
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define periodic real Theta dimension.
 */
struct Theta
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;
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
int constexpr BSDegreeP = 3;

bool constexpr BsplineOnUniformCellsR = false;
bool constexpr BsplineOnUniformCellsP = false;

struct BSplinesR
    : std::conditional_t<
              BsplineOnUniformCellsR,
              ddc::UniformBSplines<R, BSDegreeR>,
              ddc::NonUniformBSplines<R, BSDegreeR>>
{
};
struct BSplinesTheta
    : std::conditional_t<
              BsplineOnUniformCellsP,
              ddc::UniformBSplines<Theta, BSDegreeP>,
              ddc::NonUniformBSplines<Theta, BSDegreeP>>
{
};
struct PolarBSplinesRTheta : PolarBSplines<BSplinesR, BSplinesTheta, 1>
{
};

auto constexpr SplineRBoundary = ddc::BoundCond::GREVILLE;
auto constexpr SplinePBoundary = ddc::BoundCond::PERIODIC;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplinePBoundary, SplinePBoundary>;

// --- Discrete dimensions
struct GridR : SplineInterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : SplineInterpPointsTheta::interpolation_discrete_dimension_type
{
};

// --- Operators
using SplineRThetaBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        SplineRBoundary, // boundary at r=0
        SplineRBoundary, // boundary at rmax
        SplinePBoundary,
        SplinePBoundary,
        ddc::SplineSolver::LAPACK,
        GridR,
        GridTheta>;

using SplineRThetaEvaluatorConstBound = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at r=0
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        GridR,
        GridTheta>;

using SplineRThetaEvaluatorNullBound = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule, // boundary at r=0
        ddc::NullExtrapolationRule, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        GridR,
        GridTheta>;


// --- Index definitions
using IdxR = Idx<GridR>;
using IdxTheta = Idx<GridTheta>;
using IdxRTheta = Idx<GridR, GridTheta>;

// --- IVect definitions
using IdxStepR = IdxStep<GridR>;
using IdxStepTheta = IdxStep<GridTheta>;
using IdxStepRTheta = IdxStep<GridR, GridTheta>;

// --- Domain definitions
using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;

using BSIdxRangeR = IdxRange<BSplinesR>;
using BSIdxRangeTheta = IdxRange<BSplinesTheta>;
using BSIdxRangeRTheta = IdxRange<BSplinesR, BSplinesTheta>;
using BSIdxRangePolar = IdxRange<PolarBSplinesRTheta>;


// --- Chunk definitions
template <class ElementType>
using FieldMemR = host_t<FieldMem<ElementType, IdxRangeR>>;

template <class ElementType>
using FieldMemTheta = host_t<FieldMem<ElementType, IdxRangeTheta>>;

template <class ElementType>
using FieldMemRTheta = host_t<FieldMem<ElementType, IdxRangeRTheta>>;

using DFieldMemR = FieldMemR<double>;
using DFieldMemTheta = FieldMemTheta<double>;
using DFieldMemRTheta = FieldMemRTheta<double>;

// --- Span definitions
template <class ElementType>
using FieldR = host_t<Field<ElementType, IdxRangeR>>;

template <class ElementType>
using FieldTheta = host_t<Field<ElementType, IdxRangeTheta>>;

// Equivalent to host_t<Field<ElementType, IdxRangeRTheta>> but used for type deductions
template <class ElementType>
using FieldRTheta
        = Field<ElementType,
                IdxRangeRTheta,
                std::experimental::layout_right,
                Kokkos::DefaultHostExecutionSpace::memory_space>;

using DFieldR = FieldR<double>;
using DFieldTheta = FieldTheta<double>;
using DFieldRTheta = FieldRTheta<double>;

// --- View definitions
template <class ElementType>
using ConstFieldR = host_t<ConstField<ElementType const, IdxRangeR>>;

template <class ElementType>
using ConstFieldTheta = host_t<ConstField<ElementType const, IdxRangeTheta>>;

template <class ElementType>
using ConstFieldRTheta = host_t<ConstField<ElementType const, IdxRangeRTheta>>;

using DConstFieldR = ConstFieldR<double>;
using DConstFieldTheta = ConstFieldTheta<double>;
using DConstFieldRTheta = ConstFieldRTheta<double>;

// --- Spline representation definitions
using Spline2D = host_t<FieldMem<double, BSIdxRangeRTheta>>;
using Spline2DField = host_t<Field<double, BSIdxRangeRTheta>>;
using Spline2DConstField = host_t<Field<double const, BSIdxRangeRTheta>>;

/**
 * @brief Tag the polar B-splines decomposition of a function.
 *
 * Store the polar B-splines coefficients of the function.
 */
using SplinePolar = PolarSpline<PolarBSplinesRTheta>;

/**
 * @brief Type of the index of an element of polar B-splines.
 */
using IdxPolarBspl = Idx<PolarBSplinesRTheta>;


// --- VectorField definitions
template <class Dim1, class Dim2>
using DVectorFieldMemRTheta = VectorField<double, IdxRangeRTheta, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using DVectorFieldRTheta = VectorFieldSpan<double, IdxRangeRTheta, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using DConstVectorFieldRTheta = VectorFieldView<double, IdxRangeRTheta, NDTag<Dim1, Dim2>>;



template <class Dim1, class Dim2>
using VectorSplineCoeffsMem2D = VectorField<double, BSIdxRangeRTheta, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using VectorSplineCoeffs2D = VectorFieldSpan<double, BSIdxRangeRTheta, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using ConstVectorSplineCoeffs2D = VectorFieldView<double, BSIdxRangeRTheta, NDTag<Dim1, Dim2>>;



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
