// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/spline_boundary_conditions.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include <species_info.hpp>

struct RDimX
{
    static bool constexpr PERIODIC = true;
};

struct RDimY
{
    static bool constexpr PERIODIC = true;
};

struct RDimVx
{
    static bool constexpr PERIODIC = false;
};

struct RDimVy
{
    static bool constexpr PERIODIC = false;
};

using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;
using CoordXY = ddc::Coordinate<RDimX, RDimY>;

using CoordVx = ddc::Coordinate<RDimVx>;
using CoordVy = ddc::Coordinate<RDimVy>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeY = 3;

int constexpr BSDegreeVx = 3;
int constexpr BSDegreeVy = 3;

using BSplinesX = UniformBSplines<RDimX, BSDegreeX>;
using BSplinesY = UniformBSplines<RDimY, BSDegreeY>;

using BSplinesVx = UniformBSplines<RDimVx, BSDegreeVx>;
using BSplinesVy = UniformBSplines<RDimVy, BSDegreeVy>;

BoundCond constexpr SplineXBoundary = BoundCond::PERIODIC;
BoundCond constexpr SplineYBoundary = BoundCond::PERIODIC;

// IDim definition
using IDimSp = SpeciesInformation;
using InterpPointsX = GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using IDimX = typename InterpPointsX::interpolation_mesh_type;

using InterpPointsY = GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
using IDimY = typename InterpPointsY::interpolation_mesh_type;

using InterpPointsVx
        = GrevilleInterpolationPoints<BSplinesVx, BoundCond::HERMITE, BoundCond::HERMITE>;
using IDimVx = typename InterpPointsVx::interpolation_mesh_type;

using InterpPointsVy
        = GrevilleInterpolationPoints<BSplinesVy, BoundCond::HERMITE, BoundCond::HERMITE>;
using IDimVy = typename InterpPointsVy::interpolation_mesh_type;

// SplineBuilder and SplineEvaluator definition
using SplineXBuilder = SplineBuilder<BSplinesX, IDimX, SplineXBoundary, SplineXBoundary>;
using SplineYBuilder = SplineBuilder<BSplinesY, IDimY, SplineYBoundary, SplineYBoundary>;
using SplineXYBuilder = SplineBuilder2D<SplineXBuilder, SplineYBuilder>;
using SplineVxBuilder = SplineBuilder<BSplinesVx, IDimVx, BoundCond::HERMITE, BoundCond::HERMITE>;
using SplineVyBuilder = SplineBuilder<BSplinesVy, IDimVy, BoundCond::HERMITE, BoundCond::HERMITE>;
using SplineVxVyBuilder = SplineBuilder2D<SplineVxBuilder, SplineVyBuilder>;
using SplineXYEvaluator = SplineEvaluator2D<BSplinesX, BSplinesY>;
using SplineVxVyEvaluator = SplineEvaluator2D<BSplinesVx, BSplinesVy>;

using BSDomainX = ddc::DiscreteDomain<BSplinesX>;
using BSDomainY = ddc::DiscreteDomain<BSplinesY>;
using BSDomainXY = ddc::DiscreteDomain<BSplinesX, BSplinesY>;
using BSDomainVx = ddc::DiscreteDomain<BSplinesVx>;
using BSDomainVy = ddc::DiscreteDomain<BSplinesVy>;
using BSDomainVxVy = ddc::DiscreteDomain<BSplinesVx, BSplinesVy>;

template <class ElementType>
using BSViewXY = ddc::ChunkSpan<ElementType const, BSDomainXY>;
using DBSViewXY = BSViewXY<double>;

// Index
using IndexSp = ddc::DiscreteElement<IDimSp>;
using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using IndexVx = ddc::DiscreteElement<IDimVx>;
using IndexVy = ddc::DiscreteElement<IDimVy>;
using IndexXYVxVy = ddc::DiscreteElement<IDimX, IDimY, IDimVx, IDimVy>;

// IVect definition
using IVectSp = ddc::DiscreteVector<IDimSp>;
using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;
using IVectVx = ddc::DiscreteVector<IDimVx>;
using IVectVy = ddc::DiscreteVector<IDimVy>;

// Idomain definition
using IDomainSp = ddc::DiscreteDomain<IDimSp>;
using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainVx = ddc::DiscreteDomain<IDimVx>;
using IDomainVy = ddc::DiscreteDomain<IDimVy>;
using IDomainXYVxVy = ddc::DiscreteDomain<IDimX, IDimY, IDimVx, IDimVy>;
using IDomainVxVy = ddc::DiscreteDomain<IDimVx, IDimVy>;
using IDomainSpVxVy = ddc::DiscreteDomain<IDimSp, IDimVx, IDimVy>;
using IDomainSpXYVxVy = ddc::DiscreteDomain<IDimSp, IDimX, IDimY, IDimVx, IDimVy>;

// Field definition
template <class ElementType>
using FieldSp = ddc::Chunk<ElementType, IDomainSp>;
using DFieldSp = FieldSp<double>;

template <class ElementType>
using FieldX = ddc::Chunk<ElementType, IDomainX>;
using DFieldX = FieldX<double>;

template <class ElementType>
using FieldY = ddc::Chunk<ElementType, IDomainY>;
using DFieldY = FieldY<double>;

template <class ElementType>
using FieldXY = ddc::Chunk<ElementType, IDomainXY>;
using DFieldXY = FieldXY<double>;

template <class ElementType>
using FieldVx = ddc::Chunk<ElementType, IDomainVx>;

template <class ElementType>
using FieldVy = ddc::Chunk<ElementType, IDomainVy>;

template <class ElementType>
using FieldVxVy = ddc::Chunk<ElementType, IDomainVxVy>;
using DFieldVxVy = FieldVxVy<double>;

template <class ElementType>
using FieldSpVxVy = ddc::Chunk<ElementType, IDomainSpVxVy>;
using DFieldSpVxVy = FieldSpVxVy<double>;

template <class ElementType>
using FieldSpXYVxVy = ddc::Chunk<ElementType, IDomainSpXYVxVy>;
using DFieldSpXYVxVy = FieldSpXYVxVy<double>;

//  Span definition
template <class ElementType>
using SpanX = ddc::ChunkSpan<ElementType, IDomainX>;
using DSpanX = SpanX<double>;

template <class ElementType>
using SpanY = ddc::ChunkSpan<ElementType, IDomainY>;
using DSpanY = SpanY<double>;

template <class ElementType>
using SpanXY = ddc::ChunkSpan<ElementType, IDomainXY>;
using DSpanXY = SpanXY<double>;

template <class ElementType>
using SpanVx = ddc::ChunkSpan<ElementType, IDomainVx>;
using DSpanVx = SpanVx<double>;

template <class ElementType>
using SpanVy = ddc::ChunkSpan<ElementType, IDomainVy>;
using DSpanVy = SpanVy<double>;

template <class ElementType>
using SpanVxVy = ddc::ChunkSpan<ElementType, IDomainVxVy>;
using DSpanVxVy = SpanVxVy<double>;

template <class ElementType>
using SpanSpVxVy = ddc::ChunkSpan<ElementType, IDomainSpVxVy>;
using DSpanSpVxVy = SpanSpVxVy<double>;

template <class ElementType>
using SpanSpXYVxVy = ddc::ChunkSpan<ElementType, IDomainSpXYVxVy>;
using DSpanSpXYVxVy = SpanSpXYVxVy<double>;

// View definition
template <class ElementType>
using ViewSp = ddc::ChunkSpan<ElementType const, IDomainSp>;
using DViewSp = ViewSp<double>;

template <class ElementType>
using ViewX = ddc::ChunkSpan<ElementType const, IDomainX>;

template <class ElementType>
using ViewY = ddc::ChunkSpan<ElementType const, IDomainY>;

template <class ElementType>
using ViewXY = ddc::ChunkSpan<ElementType const, IDomainXY>;
using DViewXY = ViewXY<double>;

template <class ElementType>
using ViewVx = ddc::ChunkSpan<ElementType const, IDomainVx>;

template <class ElementType>
using ViewVy = ddc::ChunkSpan<ElementType const, IDomainVy>;

template <class ElementType>
using ViewSpVxVy = ddc::ChunkSpan<ElementType const, IDomainSpVxVy>;
using DViewSpVxVy = ViewSpVxVy<double>;

template <class ElementType>
using ViewSpXYVxVy = ddc::ChunkSpan<ElementType const, IDomainSpXYVxVy>;
using DViewSpXYVxVy = ViewSpXYVxVy<double>;

// For Fourier
using RDimFx = ddc::Fourier<RDimX>;
using RDimFy = ddc::Fourier<RDimY>;
using CoordFx = ddc::Coordinate<RDimFx>;
using CoordFy = ddc::Coordinate<RDimFy>;
using IDimFx = ddc::PeriodicSampling<RDimFx>;
using IDimFy = ddc::PeriodicSampling<RDimFy>;
using IndexFx = ddc::DiscreteElement<IDimFx>;
using IndexFy = ddc::DiscreteElement<IDimFy>;
using IndexFxFy = ddc::DiscreteElement<IDimFx, IDimFy>;
using IVectFx = ddc::DiscreteVector<IDimFx>;
using IVectFy = ddc::DiscreteVector<IDimFy>;
using IVectFxFy = ddc::DiscreteVector<IDimFx, IDimFy>;
using IDomainFx = ddc::DiscreteDomain<IDimFx>;
using IDomainFy = ddc::DiscreteDomain<IDimFy>;
using IDomainFxFy = ddc::DiscreteDomain<IDimFx, IDimFy>;

class GeometryXYVxVy
{
public:
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, IDimX>,
            IDimVx,
            std::conditional_t<std::is_same_v<T, IDimY>, IDimVy, void>>;

    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, IDimVx>, IDimX, std::conditional_t<std::is_same_v<T, IDimVy>, IDimY, void>>;

    using DDimSp = IDimSp;

    using SpatialDDom = IDomainXY;

    using VelocityDDom = IDomainVxVy;

    using FdistribuDDom = IDomainSpXYVxVy;
};
