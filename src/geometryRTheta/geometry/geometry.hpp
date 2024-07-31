#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/polar_bsplines.hpp>

#include <directional_tag.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>


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
struct RDimR
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
struct RDimP
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
struct RDimVr
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
struct RDimVp
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;
};


using CoordR = ddc::Coordinate<RDimR>;
using CoordP = ddc::Coordinate<RDimP>;
using CoordRP = ddc::Coordinate<RDimR, RDimP>;

using CoordVr = ddc::Coordinate<RDimVr>;
using CoordVp = ddc::Coordinate<RDimVp>;

// --- Spline definitions
int constexpr BSDegreeR = 3;
int constexpr BSDegreeP = 3;

bool constexpr BsplineOnUniformCellsR = false;
bool constexpr BsplineOnUniformCellsP = false;

struct BSplinesR
    : std::conditional_t<
              BsplineOnUniformCellsR,
              ddc::UniformBSplines<RDimR, BSDegreeR>,
              ddc::NonUniformBSplines<RDimR, BSDegreeR>>
{
};
struct BSplinesP
    : std::conditional_t<
              BsplineOnUniformCellsP,
              ddc::UniformBSplines<RDimP, BSDegreeP>,
              ddc::NonUniformBSplines<RDimP, BSDegreeP>>
{
};
struct PolarBSplinesRP : PolarBSplines<BSplinesR, BSplinesP, 1>
{
};

auto constexpr SplineRBoundary = ddc::BoundCond::GREVILLE;
auto constexpr SplinePBoundary = ddc::BoundCond::PERIODIC;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsP
        = ddc::GrevilleInterpolationPoints<BSplinesP, SplinePBoundary, SplinePBoundary>;

// --- Discrete dimensions
struct IDimR : SplineInterpPointsR::interpolation_discrete_dimension_type
{
};
struct IDimP : SplineInterpPointsP::interpolation_discrete_dimension_type
{
};

// --- Operators
using SplineRPBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesP,
        IDimR,
        IDimP,
        SplineRBoundary, // boundary at r=0
        SplineRBoundary, // boundary at rmax
        SplinePBoundary,
        SplinePBoundary,
        ddc::SplineSolver::LAPACK,
        IDimR,
        IDimP>;

using SplineRPEvaluatorConstBound = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesP,
        IDimR,
        IDimP,
        ddc::ConstantExtrapolationRule<RDimR, RDimP>, // boundary at r=0
        ddc::ConstantExtrapolationRule<RDimR, RDimP>, // boundary at rmax
        ddc::PeriodicExtrapolationRule<RDimP>,
        ddc::PeriodicExtrapolationRule<RDimP>,
        IDimR,
        IDimP>;

using SplineRPEvaluatorNullBound = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesP,
        IDimR,
        IDimP,
        ddc::NullExtrapolationRule, // boundary at r=0
        ddc::NullExtrapolationRule, // boundary at rmax
        ddc::PeriodicExtrapolationRule<RDimP>,
        ddc::PeriodicExtrapolationRule<RDimP>,
        IDimR,
        IDimP>;


// --- Index definitions
using IndexR = ddc::DiscreteElement<IDimR>;
using IndexP = ddc::DiscreteElement<IDimP>;
using IndexRP = ddc::DiscreteElement<IDimR, IDimP>;

// --- IVect definitions
using IVectR = ddc::DiscreteVector<IDimR>;
using IVectP = ddc::DiscreteVector<IDimP>;
using IVectRP = ddc::DiscreteVector<IDimR, IDimP>;

// --- Domain definitions
using IDomainR = ddc::DiscreteDomain<IDimR>;
using IDomainP = ddc::DiscreteDomain<IDimP>;
using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;

using BSDomainR = ddc::DiscreteDomain<BSplinesR>;
using BSDomainP = ddc::DiscreteDomain<BSplinesP>;
using BSDomainRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;
using BSDomainPolar = ddc::DiscreteDomain<PolarBSplinesRP>;


// --- Chunk definitions
template <class ElementType>
using FieldR = ddc::Chunk<ElementType, IDomainR>;

template <class ElementType>
using FieldP = ddc::Chunk<ElementType, IDomainP>;

template <class ElementType>
using FieldRP = ddc::Chunk<ElementType, IDomainRP>;

using DFieldR = FieldR<double>;
using DFieldP = FieldP<double>;
using DFieldRP = FieldRP<double>;

// --- Span definitions
template <class ElementType>
using SpanR = ddc::ChunkSpan<ElementType, IDomainR>;

template <class ElementType>
using SpanP = ddc::ChunkSpan<ElementType, IDomainP>;

template <class ElementType>
using SpanRP = ddc::ChunkSpan<ElementType, IDomainRP>;

using DSpanR = SpanR<double>;
using DSpanP = SpanP<double>;
using DSpanRP = SpanRP<double>;

// --- View definitions
template <class ElementType>
using ViewR = ddc::ChunkView<ElementType const, IDomainR>;

template <class ElementType>
using ViewP = ddc::ChunkView<ElementType const, IDomainP>;

template <class ElementType>
using ViewRP = ddc::ChunkView<ElementType const, IDomainRP>;

using DViewR = ViewR<double>;
using DViewP = ViewP<double>;
using DViewRP = ViewRP<double>;

// --- Spline representation definitions
using Spline2D = ddc::Chunk<double, BSDomainRP>;
using Spline2DSpan = ddc::ChunkSpan<double, BSDomainRP>;
using Spline2DView = ddc::ChunkSpan<double const, BSDomainRP>;

/**
 * @brief Tag the polar B-splines decomposition of a function.
 *
 * Store the polar B-splines coefficients of the function.
 */
using SplinePolar = PolarSpline<PolarBSplinesRP>;

/**
 * @brief Type of the index of an element of polar B-splines.
 */
using IndexPolarBspl = ddc::DiscreteElement<PolarBSplinesRP>;


// --- VectorField definitions
template <class Dim1, class Dim2>
using VectorDFieldRP = VectorField<double, IDomainRP, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using VectorDSpanRP = VectorFieldSpan<double, IDomainRP, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using VectorDViewRP = VectorFieldView<double, IDomainRP, NDTag<Dim1, Dim2>>;



template <class Dim1, class Dim2>
using VectorSpline2D = VectorField<double, BSDomainRP, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using VectorSpline2DSpan = VectorFieldSpan<double, BSDomainRP, NDTag<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using VectorSpline2DView = VectorFieldView<double, BSDomainRP, NDTag<Dim1, Dim2>>;



// CARTESIAN SPACE AND VELOCITY ------------------------------------------------------------------
// --- Continuous dimensions
/**
 * @brief Define non periodic real X dimension.
 */
struct RDimX
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
struct RDimY
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
struct RDimVx
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
struct RDimVy
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};


using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;
using CoordXY = ddc::Coordinate<RDimX, RDimY>;

using CoordVx = ddc::Coordinate<RDimVx>;
using CoordVy = ddc::Coordinate<RDimVy>;
