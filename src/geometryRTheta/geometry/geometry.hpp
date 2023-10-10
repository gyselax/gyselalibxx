#pragma once

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/polar_bsplines.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include <species_info.hpp>

/*
 * @file geometry.hpp
 *
 * Definition of
 *   - @f$ r@f$, @f$ \theta@f$, @f$(r, \theta)@f$ dimensions.
 *   - @f$x@f$, @f$y@f$, @f$(x, y)@f$ dimensions.
 */


/**
 * @brief Define non periodic X dimension.
 */
struct DimX
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define non periodic Y dimension.
 */
struct DimY
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define non periodic R dimension.
 */
struct DimR
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define periodic Theta dimension.
 */
struct DimP
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;
};

using CoordR = ddc::Coordinate<DimR>;
using CoordP = ddc::Coordinate<DimP>;
using CoordRP = ddc::Coordinate<DimR, DimP>;

int constexpr BSDegree = 3;

using BSplinesR = NonUniformBSplines<DimR, BSDegree>;
using BSplinesP = NonUniformBSplines<DimP, BSDegree>;
using PolarBSplinesRP = PolarBSplines<BSplinesR, BSplinesP, 1>;

using InterpPointsR
        = GrevilleInterpolationPoints<BSplinesR, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using InterpPointsP
        = GrevilleInterpolationPoints<BSplinesP, BoundCond::PERIODIC, BoundCond::PERIODIC>;


using IDimR = typename InterpPointsR::interpolation_mesh_type;
using IDimP = typename InterpPointsP::interpolation_mesh_type;

using SplineRBuilder = SplineBuilder<BSplinesR, IDimR, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using SplinePBuilder = SplineBuilder<BSplinesP, IDimP, BoundCond::PERIODIC, BoundCond::PERIODIC>;
using SplineRPBuilder = SplineBuilder2D<SplineRBuilder, SplinePBuilder>;

using SplineRPEvaluator = SplineEvaluator2D<BSplinesR, BSplinesP>;

using BSDomainR = ddc::DiscreteDomain<BSplinesR>;
using BSDomainP = ddc::DiscreteDomain<BSplinesP>;
using BSDomainRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;
using BSDomainPolar = ddc::DiscreteDomain<PolarBSplinesRP>;

using IndexR = ddc::DiscreteElement<IDimR>;
using IndexP = ddc::DiscreteElement<IDimP>;
using IndexRP = ddc::DiscreteElement<IDimR, IDimP>;

using IVectR = ddc::DiscreteVector<IDimR>;
using IVectP = ddc::DiscreteVector<IDimP>;
using IVectRP = ddc::DiscreteVector<IDimR, IDimP>;

using IDomainR = ddc::DiscreteDomain<IDimR>;
using IDomainP = ddc::DiscreteDomain<IDimP>;
using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;


template <class ElementType>
using FieldR = ddc::Chunk<ElementType, IDomainR>;

template <class ElementType>
using FieldP = ddc::Chunk<ElementType, IDomainP>;

template <class ElementType>
using FieldRP = ddc::Chunk<ElementType, IDomainRP>;


using DFieldR = FieldR<double>;
using DFieldP = FieldP<double>;
using DFieldRP = FieldRP<double>;


template <class ElementType>
using SpanR = ddc::ChunkSpan<ElementType, IDomainR>;

template <class ElementType>
using SpanP = ddc::ChunkSpan<ElementType, IDomainP>;

template <class ElementType>
using SpanRP = ddc::ChunkSpan<ElementType, IDomainRP>;


using DSpanR = SpanR<double>;
using DSpanP = SpanP<double>;
using DSpanRP = SpanRP<double>;


template <class ElementType>
using ViewR = ddc::ChunkView<ElementType const, IDomainR>;

template <class ElementType>
using ViewP = ddc::ChunkView<ElementType const, IDomainP>;

template <class ElementType>
using ViewRP = ddc::ChunkView<ElementType const, IDomainRP>;


using DViewR = ViewR<double>;
using DViewP = ViewP<double>;
using DViewRP = ViewRP<double>;

using Spline2D = ddc::Chunk<double, BSDomainRP>;
using Spline2DSpan = ddc::ChunkSpan<double, BSDomainRP>;
using Spline2DView = ddc::ChunkSpan<double const, BSDomainRP>;


// VELOCITY ----------------------------------------------------------
/**
 * @brief Define non periodic R velocity dimension.
 */
struct DimVr
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};

/**
 * @brief Define periodic Theta velocity dimension.
 */
struct DimVp
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;
};

using CoordVr = ddc::Coordinate<DimVr>;
using CoordVp = ddc::Coordinate<DimVp>;
using CoordVrVp = ddc::Coordinate<DimVr, DimVp>;

int constexpr BSDegreeV = 3;

using BSplinesVr = NonUniformBSplines<DimVr, BSDegreeV>;
using BSplinesVp = NonUniformBSplines<DimVp, BSDegreeV>;
using PolarBSplinesVrVp = PolarBSplines<BSplinesVr, BSplinesVp, 1>;


using InterpPointsVr
        = GrevilleInterpolationPoints<BSplinesVr, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using InterpPointsVp
        = GrevilleInterpolationPoints<BSplinesVp, BoundCond::PERIODIC, BoundCond::PERIODIC>;


using IDimVr = typename InterpPointsVr::interpolation_mesh_type;
using IDimVp = typename InterpPointsVp::interpolation_mesh_type;


using SplineVrBuilder = SplineBuilder<BSplinesVr, IDimVr, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using SplineVpBuilder = SplineBuilder<BSplinesVp, IDimVp, BoundCond::PERIODIC, BoundCond::PERIODIC>;
using SplineVrVpBuilder = SplineBuilder2D<SplineVrBuilder, SplineVpBuilder>;


using BSDomainVr = ddc::DiscreteDomain<BSplinesVr>;
using BSDomainVp = ddc::DiscreteDomain<BSplinesVp>;
using BSDomainVrVp = ddc::DiscreteDomain<BSplinesVr, BSplinesVp>;
using BSDomainPolar = ddc::DiscreteDomain<PolarBSplinesRP>;


using IndexVr = ddc::DiscreteElement<IDimVr>;
using IndexVp = ddc::DiscreteElement<IDimVp>;
using IndexVrVp = ddc::DiscreteElement<IDimVr, IDimVp>;


using IVectVr = ddc::DiscreteVector<IDimVr>;
using IVectVp = ddc::DiscreteVector<IDimVp>;
using IVectVrVp = ddc::DiscreteVector<IDimVr, IDimVp>;


using IDomainVr = ddc::DiscreteDomain<IDimVr>;
using IDomainVp = ddc::DiscreteDomain<IDimVp>;
using IDomainVrVp = ddc::DiscreteDomain<IDimVr, IDimVp>;


template <class ElementType>
using FieldVr = ddc::Chunk<ElementType, IDomainVr>;

template <class ElementType>
using FieldVp = ddc::Chunk<ElementType, IDomainVp>;

template <class ElementType>
using FieldVrVp = ddc::Chunk<ElementType, IDomainVrVp>;


using DFieldVr = FieldVr<double>;
using DFieldVp = FieldVp<double>;
using DFieldVrVp = FieldVrVp<double>;


template <class ElementType>
using SpanVr = ddc::ChunkSpan<ElementType, IDomainVr>;

template <class ElementType>
using SpanVp = ddc::ChunkSpan<ElementType, IDomainVp>;

template <class ElementType>
using SpanVrVp = ddc::ChunkSpan<ElementType, IDomainVrVp>;


using DSpanVr = SpanVr<double>;
using DSpanVp = SpanVp<double>;
using DSpanVrVp = SpanVrVp<double>;


template <class ElementType>
using ViewVr = ddc::ChunkView<ElementType const, IDomainVr>;

template <class ElementType>
using ViewVp = ddc::ChunkView<ElementType const, IDomainVp>;

template <class ElementType>
using ViewVrVp = ddc::ChunkView<ElementType const, IDomainVrVp>;


using DViewVr = ViewVr<double>;
using DViewVp = ViewVp<double>;
using DViewVrVp = ViewVrVp<double>;


using Spline2DV = ddc::Chunk<double, BSDomainVrVp>;
using Spline2DSpanV = ddc::ChunkSpan<double, BSDomainVrVp>;
using Spline2DViewV = ddc::ChunkSpan<double const, BSDomainVrVp>;



// SPACE AND VELOCITY ------------------------------------------------
using IndexRVr = ddc::DiscreteElement<IDimR, IDimVr>;
using IndexPVp = ddc::DiscreteElement<IDimP, IDimVp>;
using IndexRPVrVp = ddc::DiscreteElement<IDimR, IDimP, IDimVr, IDimVp>;


using IVectRVr = ddc::DiscreteVector<IDimR, IDimVr>;
using IVectPVp = ddc::DiscreteVector<IDimP, IDimVp>;
using IVectRPVrVp = ddc::DiscreteVector<IDimR, IDimP, IDimVr, IDimVp>;


using IDomainRVr = ddc::DiscreteDomain<IDimR, IDimVr>;
using IDomainPVp = ddc::DiscreteDomain<IDimP, IDimVp>;
using IDomainRPVrVp = ddc::DiscreteDomain<IDimR, IDimP, IDimVr, IDimVp>;


template <class ElementType>
using FieldRVr = ddc::Chunk<ElementType, IDomainRVr>;

template <class ElementType>
using FieldPVp = ddc::Chunk<ElementType, IDomainPVp>;

template <class ElementType>
using FieldRPVrVp = ddc::Chunk<ElementType, IDomainRPVrVp>;


using DFieldRVr = FieldRVr<double>;
using DFieldPVp = FieldPVp<double>;
using DFieldRPVrVp = FieldRPVrVp<double>;


template <class ElementType>
using SpanRVr = ddc::ChunkSpan<ElementType, IDomainRVr>;

template <class ElementType>
using SpanPVp = ddc::ChunkSpan<ElementType, IDomainPVp>;

template <class ElementType>
using SpanRPVrVp = ddc::ChunkSpan<ElementType, IDomainRPVrVp>;


using DSpanRVr = SpanRVr<double>;
using DSpanPVp = SpanPVp<double>;
using DSpanRPVrVp = SpanRPVrVp<double>;


template <class ElementType>
using ViewRVr = ddc::ChunkView<ElementType const, IDomainRVr>;

template <class ElementType>
using ViewPVp = ddc::ChunkView<ElementType const, IDomainPVp>;

template <class ElementType>
using ViewRPVrVp = ddc::ChunkView<ElementType const, IDomainRPVrVp>;


using DViewRVr = ViewRVr<double>;
using DViewPVp = ViewPVp<double>;
using DViewRPVrVp = ViewRPVrVp<double>;



// CARTESIAN SPACE AND VELOCITY -------------------------------------
// -- X DIMENSION ---------------------------------------------------
/**
 * @brief Define non periodic X velocity dimension.
 */
struct DimVx
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};



using CoordX = ddc::Coordinate<DimX>;
using CoordVx = ddc::Coordinate<DimVx>;
using CoordXVx = ddc::Coordinate<DimX, DimVx>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

using BSplinesX = UniformBSplines<DimX, BSDegreeX>;
using BSplinesVx = UniformBSplines<DimVx, BSDegreeVx>;


auto constexpr SplineXBoundary = DimX::PERIODIC ? BoundCond::PERIODIC : BoundCond::GREVILLE;
auto constexpr SplineVxBoundary = BoundCond::HERMITE;

using InterpPointsX = GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using IDimX = typename InterpPointsX::interpolation_mesh_type;
using SplineXBuilder = SplineBuilder<BSplinesX, IDimX, SplineXBoundary, SplineXBoundary>;

using InterpPointsVx = GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;
using IDimVx = typename InterpPointsVx::interpolation_mesh_type;
using SplineVxBuilder = SplineBuilder<BSplinesVx, IDimVx, SplineVxBoundary, SplineVxBoundary>;

// Species dimension
using IDimSp = SpeciesInformation;


using IndexX = ddc::DiscreteElement<IDimX>;
using IndexVx = ddc::DiscreteElement<IDimVx>;
using IndexXVx = ddc::DiscreteElement<IDimX, IDimVx>;

using IndexSp = ddc::DiscreteElement<IDimSp>;
using IndexSpX = ddc::DiscreteElement<IDimSp, IDimX>;
using IndexSpVx = ddc::DiscreteElement<IDimSp, IDimVx>;
using IndexSpXVx = ddc::DiscreteElement<IDimSp, IDimX, IDimVx>;



using IVectX = ddc::DiscreteVector<IDimX>;
using IVectVx = ddc::DiscreteVector<IDimVx>;
using IVectXVx = ddc::DiscreteVector<IDimX, IDimVx>;

using IVectSp = ddc::DiscreteVector<IDimSp>;
using IVectSpX = ddc::DiscreteVector<IDimSp, IDimX>;
using IVectSpVx = ddc::DiscreteVector<IDimSp, IDimVx>;
using IVectSpXVx = ddc::DiscreteVector<IDimSp, IDimX, IDimVx>;



using BSDomainX = ddc::DiscreteDomain<BSplinesX>;
using BSDomainVx = ddc::DiscreteDomain<BSplinesVx>;

using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainVx = ddc::DiscreteDomain<IDimVx>;
using IDomainXVx = ddc::DiscreteDomain<IDimX, IDimVx>;

using IDomainSp = ddc::DiscreteDomain<IDimSp>;
using IDomainSpX = ddc::DiscreteDomain<IDimSp, IDimX>;
using IDomainSpVx = ddc::DiscreteDomain<IDimSp, IDimVx>;
using IDomainSpXVx = ddc::DiscreteDomain<IDimSp, IDimX, IDimVx>;



template <class ElementType>
using FieldX = ddc::Chunk<ElementType, IDomainX>;

template <class ElementType>
using FieldVx = ddc::Chunk<ElementType, IDomainVx>;

template <class ElementType>
using FieldSp = ddc::Chunk<ElementType, IDomainSp>;

template <class ElementType>
using FieldSpX = ddc::Chunk<ElementType, IDomainSpX>;

template <class ElementType>
using FieldSpVx = ddc::Chunk<ElementType, IDomainSpVx>;

template <class ElementType>
using FieldSpXVx = ddc::Chunk<ElementType, IDomainSpXVx>;



using DFieldX = FieldX<double>;
using DFieldVx = FieldVx<double>;

using DFieldSp = FieldSp<double>;
using DFieldSpX = FieldSpX<double>;
using DFieldSpVx = FieldSpVx<double>;
using DFieldSpXVx = FieldSpXVx<double>;



template <class ElementType>
using SpanX = ddc::ChunkSpan<ElementType, IDomainX>;

template <class ElementType>
using SpanVx = ddc::ChunkSpan<ElementType, IDomainVx>;

template <class ElementType>
using SpanSp = ddc::ChunkSpan<ElementType, IDomainSp>;

template <class ElementType>
using SpanSpX = ddc::ChunkSpan<ElementType, IDomainSpX>;

template <class ElementType>
using SpanSpVx = ddc::ChunkSpan<ElementType, IDomainSpVx>;

template <class ElementType>
using SpanSpXVx = ddc::ChunkSpan<ElementType, IDomainSpXVx>;



using DSpanX = SpanX<double>;
using DSpanVx = SpanVx<double>;

using DSpanSp = SpanSp<double>;
using DSpanSpX = SpanSpX<double>;
using DSpanSpVx = SpanSpVx<double>;
using DSpanSpXVx = SpanSpXVx<double>;



template <class ElementType>
using ViewX = ddc::ChunkSpan<ElementType const, IDomainX>;

template <class ElementType>
using ViewVx = ddc::ChunkSpan<ElementType const, IDomainVx>;

template <class ElementType>
using ViewSp = ddc::ChunkSpan<ElementType const, IDomainSp>;

template <class ElementType>
using ViewSpX = ddc::ChunkSpan<ElementType const, IDomainSpX>;

template <class ElementType>
using ViewSpVx = ddc::ChunkSpan<ElementType const, IDomainSpVx>;

template <class ElementType>
using ViewSpXVx = ddc::ChunkSpan<ElementType const, IDomainSpXVx>;

template <class ElementType>
using BSViewX = ddc::ChunkSpan<ElementType const, BSDomainX>;



using DViewX = ViewX<double>;
using DViewVx = ViewVx<double>;

using DViewSp = ViewSp<double>;
using DViewSpX = ViewSpX<double>;
using DViewSpVx = ViewSpVx<double>;
using DViewSpXVx = ViewSpXVx<double>;

using DBSViewX = BSViewX<double>;


// -- Y DIMENSION ---------------------------------------------------
/**
 * @brief Define non periodic Y velocity dimension.
 */
struct DimVy
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};


using CoordY = ddc::Coordinate<DimY>;
using CoordVy = ddc::Coordinate<DimVy>;
using CoordYVy = ddc::Coordinate<DimY, DimVy>;

int constexpr BSDegreeY = 3;
int constexpr BSDegreeVy = 3;

using BSplinesY = UniformBSplines<DimY, BSDegreeY>;
using BSplinesVy = UniformBSplines<DimVy, BSDegreeVy>;


auto constexpr SplineYBoundary = DimY::PERIODIC ? BoundCond::PERIODIC : BoundCond::GREVILLE;
auto constexpr SplineVyBoundary = BoundCond::HERMITE;

using InterpPointsY = GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
using IDimY = typename InterpPointsY::interpolation_mesh_type;
using SplineYBuilder = SplineBuilder<BSplinesY, IDimY, SplineYBoundary, SplineYBoundary>;

using InterpPointsVy = GrevilleInterpolationPoints<BSplinesVy, SplineVyBoundary, SplineVyBoundary>;
using IDimVy = typename InterpPointsVy::interpolation_mesh_type;
using SplineVyBuilder = SplineBuilder<BSplinesVy, IDimVy, SplineVyBoundary, SplineVyBoundary>;

// Species dimension
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexVy = ddc::DiscreteElement<IDimVy>;
using IndexYVy = ddc::DiscreteElement<IDimY, IDimVy>;

using IndexSpY = ddc::DiscreteElement<IDimSp, IDimY>;
using IndexSpVy = ddc::DiscreteElement<IDimSp, IDimVy>;
using IndexSpYVy = ddc::DiscreteElement<IDimSp, IDimY, IDimVy>;



using IVectY = ddc::DiscreteVector<IDimY>;
using IVectVy = ddc::DiscreteVector<IDimVy>;
using IVectYVy = ddc::DiscreteVector<IDimY, IDimVy>;

using IVectSpY = ddc::DiscreteVector<IDimSp, IDimY>;
using IVectSpVy = ddc::DiscreteVector<IDimSp, IDimVy>;
using IVectSpYVy = ddc::DiscreteVector<IDimSp, IDimY, IDimVy>;



using BSDomainY = ddc::DiscreteDomain<BSplinesY>;
using BSDomainVy = ddc::DiscreteDomain<BSplinesVy>;

using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainVy = ddc::DiscreteDomain<IDimVy>;
using IDomainYVy = ddc::DiscreteDomain<IDimY, IDimVy>;

using IDomainSpY = ddc::DiscreteDomain<IDimSp, IDimY>;
using IDomainSpVy = ddc::DiscreteDomain<IDimSp, IDimVy>;
using IDomainSpYVy = ddc::DiscreteDomain<IDimSp, IDimY, IDimVy>;



template <class ElementType>
using FieldY = ddc::Chunk<ElementType, IDomainY>;

template <class ElementType>
using FieldVy = ddc::Chunk<ElementType, IDomainVy>;

template <class ElementType>
using FieldSpY = ddc::Chunk<ElementType, IDomainSpY>;

template <class ElementType>
using FieldSpVy = ddc::Chunk<ElementType, IDomainSpVy>;

template <class ElementType>
using FieldSpYVy = ddc::Chunk<ElementType, IDomainSpYVy>;



using DFieldY = FieldY<double>;
using DFieldVy = FieldVy<double>;

using DFieldSpY = FieldSpY<double>;
using DFieldSpVy = FieldSpVy<double>;
using DFieldSpYVy = FieldSpYVy<double>;



template <class ElementType>
using SpanY = ddc::ChunkSpan<ElementType, IDomainY>;

template <class ElementType>
using SpanVy = ddc::ChunkSpan<ElementType, IDomainVy>;

template <class ElementType>
using SpanSpY = ddc::ChunkSpan<ElementType, IDomainSpY>;

template <class ElementType>
using SpanSpVy = ddc::ChunkSpan<ElementType, IDomainSpVy>;

template <class ElementType>
using SpanSpYVy = ddc::ChunkSpan<ElementType, IDomainSpYVy>;



using DSpanY = SpanY<double>;
using DSpanVy = SpanVy<double>;

using DSpanSpY = SpanSpY<double>;
using DSpanSpVy = SpanSpVy<double>;
using DSpanSpYVy = SpanSpYVy<double>;



template <class ElementType>
using ViewY = ddc::ChunkSpan<ElementType const, IDomainY>;

template <class ElementType>
using ViewVy = ddc::ChunkSpan<ElementType const, IDomainVy>;

template <class ElementType>
using ViewSpY = ddc::ChunkSpan<ElementType const, IDomainSpY>;

template <class ElementType>
using ViewSpVy = ddc::ChunkSpan<ElementType const, IDomainSpVy>;

template <class ElementType>
using ViewSpYVy = ddc::ChunkSpan<ElementType const, IDomainSpYVy>;

template <class ElementType>
using BSViewY = ddc::ChunkSpan<ElementType const, BSDomainY>;



using DViewY = ViewY<double>;
using DViewVy = ViewVy<double>;

using DViewSpY = ViewSpY<double>;
using DViewSpVy = ViewSpVy<double>;
using DViewSpYVy = ViewSpYVy<double>;

using DBSViewY = BSViewY<double>;



// -- X Y DIMENSIONS ------------------------------------------------
using CoordXY = ddc::Coordinate<DimX, DimY>;
using CoordVxVy = ddc::Coordinate<DimVx, DimVy>;
using CoordXYVxVy = ddc::Coordinate<DimX, DimY, DimVx, DimVy>;


using SplineXYBuilder = SplineBuilder2D<SplineXBuilder, SplineYBuilder>;
using SplineVxVyBuilder = SplineBuilder2D<SplineVxBuilder, SplineVyBuilder>;


// Species dimension
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using IndexVxVy = ddc::DiscreteElement<IDimVx, IDimVy>;
using IndexXYVxVy = ddc::DiscreteElement<IDimX, IDimY, IDimVx, IDimVy>;

using IndexSpXY = ddc::DiscreteElement<IDimSp, IDimX, IDimY>;
using IndexSpVxVy = ddc::DiscreteElement<IDimSp, IDimVx, IDimVy>;
using IndexSpXYVxVy = ddc::DiscreteElement<IDimSp, IDimX, IDimY, IDimVx, IDimVy>;



using IVectXY = ddc::DiscreteVector<IDimX, IDimY>;
using IVectVxVy = ddc::DiscreteVector<IDimVx, IDimVy>;
using IVectXYVxVy = ddc::DiscreteVector<IDimX, IDimY, IDimVx, IDimVy>;

using IVectSpXY = ddc::DiscreteVector<IDimSp, IDimX, IDimY>;
using IVectSpVxVy = ddc::DiscreteVector<IDimSp, IDimVx, IDimVy>;
using IVectSpXYVxVy = ddc::DiscreteVector<IDimSp, IDimX, IDimY, IDimVx, IDimVy>;



using BSDomainXY = ddc::DiscreteDomain<BSplinesX, BSplinesY>;
using BSDomainVxVy = ddc::DiscreteDomain<BSplinesVx, BSplinesVy>;

using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainVxVy = ddc::DiscreteDomain<IDimVx, IDimVy>;
using IDomainXYVxVy = ddc::DiscreteDomain<IDimX, IDimY, IDimVx, IDimVy>;

using IDomainSpXY = ddc::DiscreteDomain<IDimSp, IDimX, IDimY>;
using IDomainSpVxVy = ddc::DiscreteDomain<IDimSp, IDimVx, IDimVy>;
using IDomainSpXYVxVy = ddc::DiscreteDomain<IDimSp, IDimX, IDimY, IDimVx, IDimVy>;



template <class ElementType>
using FieldXY = ddc::Chunk<ElementType, IDomainXY>;

template <class ElementType>
using FieldVxVy = ddc::Chunk<ElementType, IDomainVxVy>;

template <class ElementType>
using FieldSpXY = ddc::Chunk<ElementType, IDomainSpXY>;

template <class ElementType>
using FieldSpVxVy = ddc::Chunk<ElementType, IDomainSpVxVy>;

template <class ElementType>
using FieldSpXYVxVy = ddc::Chunk<ElementType, IDomainSpXYVxVy>;



using DFieldXY = FieldXY<double>;
using DFieldVxVy = FieldVxVy<double>;

using DFieldSpVxVy = FieldSpVxVy<double>;
using DFieldSpXYVxVy = FieldSpXYVxVy<double>;



template <class ElementType>
using SpanXY = ddc::ChunkSpan<ElementType, IDomainXY>;

template <class ElementType>
using SpanVxVy = ddc::ChunkSpan<ElementType, IDomainVy>;

template <class ElementType>
using SpanSpXY = ddc::ChunkSpan<ElementType, IDomainSpXY>;

template <class ElementType>
using SpanSpVxVy = ddc::ChunkSpan<ElementType, IDomainSpVy>;

template <class ElementType>
using SpanSpXYVxVy = ddc::ChunkSpan<ElementType, IDomainSpXYVxVy>;



using DSpanXY = SpanXY<double>;
using DSpanVxVy = SpanVxVy<double>;

using DSpanSpXY = SpanSpXY<double>;
using DSpanSpVxVy = SpanSpVxVy<double>;
using DSpanSpXYVxVy = SpanSpXYVxVy<double>;



template <class ElementType>
using ViewXY = ddc::ChunkSpan<ElementType const, IDomainXY>;

template <class ElementType>
using ViewVxVy = ddc::ChunkSpan<ElementType const, IDomainVxVy>;

template <class ElementType>
using ViewSpXY = ddc::ChunkSpan<ElementType const, IDomainSpXY>;

template <class ElementType>
using ViewSpVxVy = ddc::ChunkSpan<ElementType const, IDomainSpVxVy>;

template <class ElementType>
using ViewSpXYVxVy = ddc::ChunkSpan<ElementType const, IDomainSpXYVxVy>;

template <class ElementType>
using BSViewXY = ddc::ChunkSpan<ElementType const, BSDomainXY>;



using DViewXY = ViewXY<double>;
using DViewVxVy = ViewVxVy<double>;

using DViewSpXY = ViewSpXY<double>;
using DViewSpVxVy = ViewSpVxVy<double>;
using DViewSpXYVxVy = ViewSpXYVxVy<double>;

using DBSViewXY = BSViewXY<double>;
