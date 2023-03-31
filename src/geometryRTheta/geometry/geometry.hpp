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

using BSDomainR = ddc::DiscreteDomain<BSplinesR>;
using BSDomainR = ddc::DiscreteDomain<BSplinesR>;
using BSDomainP = ddc::DiscreteDomain<BSplinesP>;
using BSDomainP = ddc::DiscreteDomain<BSplinesP>;
using BSDomainRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;
using BSDomainRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;
using BSDomainPolar = ddc::DiscreteDomain<PolarBSplinesRP>;
using BSDomainPolar = ddc::DiscreteDomain<PolarBSplinesRP>;

using IndexR = ddc::DiscreteElement<IDimR>;
using IndexP = ddc::DiscreteElement<IDimP>;
using IndexRP = ddc::DiscreteElement<IDimR, IDimP>;

using IVectR = ddc::DiscreteVector<IDimR>;
using IVectP = ddc::DiscreteVector<IDimP>;
using IVectRP = ddc::DiscreteVector<IDimR, IDimP>;

using IDomainR = ddc::DiscreteDomain<IDimR>;
using IDomainR = ddc::DiscreteDomain<IDimR>;
using IDomainP = ddc::DiscreteDomain<IDimP>;
using IDomainP = ddc::DiscreteDomain<IDimP>;
using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;
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
