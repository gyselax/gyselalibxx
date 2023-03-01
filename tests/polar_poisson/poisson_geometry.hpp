#pragma once

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/polar_bsplines.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator.hpp>
#include <sll/spline_evaluator_2d.hpp>

struct DimX
{
    static bool constexpr PERIODIC = false;
};
struct DimY
{
    static bool constexpr PERIODIC = false;
};
struct DimR
{
    static bool constexpr PERIODIC = false;
};

struct DimP
{
    static bool constexpr PERIODIC = true;
};

using CoordR = Coordinate<DimR>;
using CoordP = Coordinate<DimP>;
using CoordRP = Coordinate<DimR, DimP>;

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

using BSDomainR = DiscreteDomain<BSplinesR>;
using BSDomainP = DiscreteDomain<BSplinesP>;
using BSDomainRP = DiscreteDomain<BSplinesR, BSplinesP>;
using BSDomainPolar = DiscreteDomain<PolarBSplinesRP>;

using IndexR = DiscreteElement<IDimR>;
using IndexP = DiscreteElement<IDimP>;
using IndexRP = DiscreteElement<IDimR, IDimP>;

using IVectR = DiscreteVector<IDimR>;
using IVectP = DiscreteVector<IDimP>;
using IVectRP = DiscreteVector<IDimR, IDimP>;

using IDomainR = DiscreteDomain<IDimR>;
using IDomainP = DiscreteDomain<IDimP>;
using IDomainRP = DiscreteDomain<IDimR, IDimP>;

template <class ElementType>
using FieldR = Chunk<ElementType, IDomainR>;

template <class ElementType>
using FieldP = Chunk<ElementType, IDomainP>;

template <class ElementType>
using FieldRP = Chunk<ElementType, IDomainRP>;

using DFieldR = FieldR<double>;
using DFieldP = FieldP<double>;
using DFieldRP = FieldRP<double>;

template <class ElementType>
using SpanR = ChunkSpan<ElementType, IDomainR>;

template <class ElementType>
using SpanP = ChunkSpan<ElementType, IDomainP>;

template <class ElementType>
using SpanRP = ChunkSpan<ElementType, IDomainRP>;

using DSpanR = SpanR<double>;
using DSpanP = SpanP<double>;
using DSpanRP = SpanRP<double>;

template <class ElementType>
using ViewR = ChunkView<ElementType const, IDomainR>;

template <class ElementType>
using ViewP = ChunkView<ElementType const, IDomainP>;

template <class ElementType>
using ViewRP = ChunkView<ElementType const, IDomainRP>;

using DViewR = ViewR<double>;
using DViewP = ViewP<double>;
using DViewRP = ViewRP<double>;

using Spline2D = Chunk<double, BSDomainRP>;
using Spline2DSpan = ChunkSpan<double, BSDomainRP>;
using Spline2DView = ChunkSpan<double const, BSDomainRP>;
