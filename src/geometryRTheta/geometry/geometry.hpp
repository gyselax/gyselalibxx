#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/polar_bsplines.hpp>

#include <directional_tag.hpp>
#include <species_info.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>


/*
 * @file geometry.hpp
 *
 * Definition of
 *   - @f$ r@f$, @f$ \theta@f$, @f$(r, \theta)@f$ dimensions.
 *   - @f$x@f$, @f$y@f$, @f$(x, y)@f$ dimensions.
 */



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

using CoordR = ddc::Coordinate<RDimR>;
using CoordP = ddc::Coordinate<RDimP>;
using CoordRP = ddc::Coordinate<RDimR, RDimP>;

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

bool constexpr UniformMeshR = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsR,
        SplineRBoundary,
        SplineRBoundary,
        BSDegreeR);
bool constexpr UniformMeshP = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsP,
        SplinePBoundary,
        SplinePBoundary,
        BSDegreeP);

struct IDimR
    : std::conditional_t<
              UniformMeshR,
              ddc::UniformPointSampling<RDimR>,
              ddc::NonUniformPointSampling<RDimR>>
{
};
struct IDimP
    : std::conditional_t<
              UniformMeshP,
              ddc::UniformPointSampling<RDimP>,
              ddc::NonUniformPointSampling<RDimP>>
{
};

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsP
        = ddc::GrevilleInterpolationPoints<BSplinesP, SplinePBoundary, SplinePBoundary>;

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
        ddc::SplineSolver::GINKGO,
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

/**
 * @brief Tag the polar B-splines decomposition of a function.
 *
 * Store the polar B-splines coefficients of the function.
 */
using SplinePolar = PolarSpline<PolarBSplinesRP>;

/**
 * @brief Tag type of an element of a couple of a B-splines in the first
 * dimension and in the second dimension.
 */
using IDimBSpline2D = ddc::DiscreteElement<BSplinesR, BSplinesP>;
/**
 * @brief Tag type of an element of polar B-splines.
 */
using IDimPolarBspl = ddc::DiscreteElement<PolarBSplinesRP>;



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


// VELOCITY ----------------------------------------------------------
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

using CoordVr = ddc::Coordinate<RDimVr>;
using CoordVp = ddc::Coordinate<RDimVp>;
using CoordVrVp = ddc::Coordinate<RDimVr, RDimVp>;

int constexpr BSDegreeVr = 3;
int constexpr BSDegreeVp = 3;

bool constexpr BsplineOnUniformCellsVr = false;
bool constexpr BsplineOnUniformCellsVp = false;

struct BSplinesVr
    : std::conditional_t<
              BsplineOnUniformCellsVr,
              ddc::UniformBSplines<RDimVr, BSDegreeVr>,
              ddc::NonUniformBSplines<RDimVr, BSDegreeVr>>
{
};
struct BSplinesVp
    : std::conditional_t<
              BsplineOnUniformCellsVp,
              ddc::UniformBSplines<RDimVp, BSDegreeVp>,
              ddc::NonUniformBSplines<RDimVp, BSDegreeVp>>
{
};
using PolarBSplinesVrVp = PolarBSplines<BSplinesVr, BSplinesVp, 1>;



auto constexpr SplineVrBoundary
        = RDimVr::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;
auto constexpr SplineVpBoundary
        = RDimVp::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;

bool constexpr UniformMeshVr = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsVr,
        SplineVrBoundary,
        SplineVrBoundary,
        BSDegreeVr);
bool constexpr UniformMeshVp = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsVp,
        SplineVpBoundary,
        SplineVpBoundary,
        BSDegreeVp);

struct IDimVr
    : std::conditional_t<
              UniformMeshVr,
              ddc::UniformPointSampling<RDimVr>,
              ddc::NonUniformPointSampling<RDimVr>>
{
};
struct IDimVp
    : std::conditional_t<
              UniformMeshVp,
              ddc::UniformPointSampling<RDimVp>,
              ddc::NonUniformPointSampling<RDimVp>>
{
};


using SplineInterpPointsVr
        = ddc::GrevilleInterpolationPoints<BSplinesVr, SplineVrBoundary, SplineVrBoundary>;
using SplineInterpPointsVp
        = ddc::GrevilleInterpolationPoints<BSplinesVp, SplineVpBoundary, SplineVpBoundary>;

using SplineVrBuilder = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVr,
        IDimVr,
        SplineVrBoundary,
        SplineVrBoundary,
        ddc::SplineSolver::GINKGO,
        IDimVr>;
using SplineVpBuilder = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVp,
        IDimVp,
        SplineVpBoundary,
        SplineVpBoundary,
        ddc::SplineSolver::GINKGO,
        IDimVp>;
using SplineVrVpBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVr,
        BSplinesVp,
        IDimVr,
        IDimVp,
        SplineVrBoundary,
        SplineVrBoundary,
        SplineVpBoundary,
        SplineVpBoundary,
        ddc::SplineSolver::GINKGO,
        IDimVr,
        IDimVp>;


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



using CoordX = ddc::Coordinate<RDimX>;
using CoordVx = ddc::Coordinate<RDimVx>;
using CoordXVx = ddc::Coordinate<RDimX, RDimVx>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

bool constexpr BsplineOnUniformCellsX = false;
bool constexpr BsplineOnUniformCellsVx = false;

struct BSplinesX
    : std::conditional_t<
              BsplineOnUniformCellsX,
              ddc::UniformBSplines<RDimX, BSDegreeX>,
              ddc::NonUniformBSplines<RDimX, BSDegreeX>>
{
};
struct BSplinesVx
    : std::conditional_t<
              BsplineOnUniformCellsVx,
              ddc::UniformBSplines<RDimVx, BSDegreeVx>,
              ddc::NonUniformBSplines<RDimVx, BSDegreeVx>>
{
};

auto constexpr SplineXBoundary
        = RDimX::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;
auto constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;


bool constexpr UniformMeshX = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsX,
        SplineXBoundary,
        SplineXBoundary,
        BSDegreeX);
bool constexpr UniformMeshVx = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsVx,
        SplineVxBoundary,
        SplineVxBoundary,
        BSDegreeVx);

struct IDimX
    : std::conditional_t<
              UniformMeshX,
              ddc::UniformPointSampling<RDimX>,
              ddc::NonUniformPointSampling<RDimX>>
{
};
struct IDimVx
    : std::conditional_t<
              UniformMeshVx,
              ddc::UniformPointSampling<RDimVx>,
              ddc::NonUniformPointSampling<RDimVx>>
{
};

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;



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


using CoordY = ddc::Coordinate<RDimY>;
using CoordVy = ddc::Coordinate<RDimVy>;
using CoordYVy = ddc::Coordinate<RDimY, RDimVy>;

int constexpr BSDegreeY = 3;
int constexpr BSDegreeVy = 3;


bool constexpr BsplineOnUniformCellsY = false;
bool constexpr BsplineOnUniformCellsVy = false;

using BSplinesY = std::conditional_t<
        BsplineOnUniformCellsY,
        ddc::UniformBSplines<RDimY, BSDegreeY>,
        ddc::NonUniformBSplines<RDimY, BSDegreeY>>;

using BSplinesVy = std::conditional_t<
        BsplineOnUniformCellsVy,
        ddc::UniformBSplines<RDimVy, BSDegreeVy>,
        ddc::NonUniformBSplines<RDimVy, BSDegreeVy>>;



auto constexpr SplineYBoundary
        = RDimY::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;
auto constexpr SplineVyBoundary = ddc::BoundCond::HERMITE;

bool constexpr UniformMeshY = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsY,
        SplineYBoundary,
        SplineYBoundary,
        BSDegreeY);
bool constexpr UniformMeshVy = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsVy,
        SplineVyBoundary,
        SplineVyBoundary,
        BSDegreeVy);

using IDimY = std::conditional_t<
        UniformMeshY,
        ddc::UniformPointSampling<RDimY>,
        ddc::NonUniformPointSampling<RDimY>>;
using IDimVy = std::conditional_t<
        UniformMeshVy,
        ddc::UniformPointSampling<RDimVy>,
        ddc::NonUniformPointSampling<RDimVy>>;

using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
using SplineInterpPointsVy
        = ddc::GrevilleInterpolationPoints<BSplinesVy, SplineVyBoundary, SplineVyBoundary>;


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
using CoordXY = ddc::Coordinate<RDimX, RDimY>;
using CoordVxVy = ddc::Coordinate<RDimVx, RDimVy>;
using CoordXYVxVy = ddc::Coordinate<RDimX, RDimY, RDimVx, RDimVy>;


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
