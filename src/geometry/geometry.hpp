// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/spline_builder.hpp>

#include <species_info.hpp>

template <class Tag>
struct Fourier
{
    using base_tag_type = Tag;
    static bool constexpr PERIODIC = Tag::PERIODIC;
};



struct RDimX
{
#ifdef PERIODIC_RDIMX
    static bool constexpr PERIODIC = true;
#else
    static bool constexpr PERIODIC = false;
#endif
};

struct RDimVx
{
    static bool constexpr PERIODIC = false;
};

struct RDimT
{
    static bool constexpr PERIODIC = false;
};



using CoordT = Coordinate<RDimT>;

using CoordX = Coordinate<RDimX>;

using CoordVx = Coordinate<RDimVx>;

using CoordXVx = Coordinate<RDimX, RDimVx>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

using BSplinesX = UniformBSplines<RDimX, BSDegreeX>;

using BSplinesVx = UniformBSplines<RDimVx, BSDegreeVx>;

auto constexpr SplineXBoundary = RDimX::PERIODIC ? BoundCond::PERIODIC : BoundCond::GREVILLE;
using SplineXBuilder = SplineBuilder<BSplinesX, SplineXBoundary, SplineXBoundary>;
using IDimX = typename SplineXBuilder::interpolation_mesh_type;

using SplineVxBuilder = SplineBuilder<BSplinesVx, BoundCond::HERMITE, BoundCond::HERMITE>;
using IDimVx = typename SplineVxBuilder::interpolation_mesh_type;

// Species dimension
using IDimSp = SpeciesInformation;


using IndexX = DiscreteElement<IDimX>;

using IndexVx = DiscreteElement<IDimVx>;

using IndexSp = DiscreteElement<IDimSp>;

using IndexSpX = DiscreteElement<IDimSp, IDimX>;

using IndexXVx = DiscreteElement<IDimX, IDimVx>;

using IndexSpVx = DiscreteElement<IDimSp, IDimVx>;

using IndexSpXVx = DiscreteElement<IDimSp, IDimX, IDimVx>;



using IVectX = DiscreteVector<IDimX>;

using IVectVx = DiscreteVector<IDimVx>;

using IVectXVx = DiscreteVector<IDimX, IDimVx>;

using IVectSpXVx = DiscreteVector<IDimSp, IDimX, IDimVx>;

using IVectSp = DiscreteVector<IDimSp>;

using IVectSpVx = DiscreteVector<IDimSp, IDimVx>;



using BSDomainX = DiscreteDomain<BSplinesX>;

using BSDomainVx = DiscreteDomain<BSplinesVx>;

using IDomainX = DiscreteDomain<IDimX>;

using IDomainVx = DiscreteDomain<IDimVx>;

using IDomainXVx = DiscreteDomain<IDimX, IDimVx>;

using IDomainSp = DiscreteDomain<IDimSp>;

using IDomainSpX = DiscreteDomain<IDimSp, IDimX>;

using IDomainSpVx = DiscreteDomain<IDimSp, IDimVx>;

using IDomainSpXVx = DiscreteDomain<IDimSp, IDimX, IDimVx>;



template <class ElementType>
using FieldX = Chunk<ElementType, IDomainX>;

template <class ElementType>
using FieldVx = Chunk<ElementType, IDomainVx>;

template <class ElementType>
using FieldSp = Chunk<ElementType, IDomainSp>;

template <class ElementType>
using FieldSpX = Chunk<ElementType, IDomainSpX>;

template <class ElementType>
using FieldSpVx = Chunk<ElementType, IDomainSpVx>;

template <class ElementType>
using FieldSpXVx = Chunk<ElementType, IDomainSpXVx>;



using DFieldX = FieldX<double>;

using DFieldVx = FieldVx<double>;

using DFieldSp = FieldSp<double>;

using DFieldSpX = FieldSpX<double>;

using DFieldSpVx = FieldSpVx<double>;

using DFieldSpXVx = FieldSpXVx<double>;


template <class ElementType>
using SpanSp = ChunkSpan<ElementType, IDomainSp>;

template <class ElementType>
using SpanX = ChunkSpan<ElementType, IDomainX>;

template <class ElementType>
using SpanVx = ChunkSpan<ElementType, IDomainVx>;

template <class ElementType>
using SpanSpXVx = ChunkSpan<ElementType, IDomainSpXVx>;

template <class ElementType>
using SpanSpVx = ChunkSpan<ElementType, IDomainSpVx>;



using DSpanSp = SpanSp<double>;

using DSpanX = SpanX<double>;

using DSpanVx = SpanVx<double>;

using DSpanSpXVx = SpanSpXVx<double>;

using DSpanSpVx = SpanSpVx<double>;



template <class ElementType>
using ViewX = ChunkSpan<ElementType const, IDomainX>;

template <class ElementType>
using ViewVx = ChunkSpan<ElementType const, IDomainVx>;

template <class ElementType>
using ViewSp = ChunkSpan<ElementType const, IDomainSp>;
template <class ElementType>
using ViewSpVx = ChunkSpan<ElementType const, IDomainSpVx>;

template <class ElementType>
using ViewSpXVx = ChunkSpan<ElementType const, IDomainSpXVx>;

template <class ElementType>
using BSViewX = ChunkSpan<ElementType const, BSDomainX>;



using DViewX = ViewX<double>;

using DViewVx = ViewVx<double>;

using DViewSp = ViewSp<double>;

using DViewSpVx = ViewSpVx<double>;

using DViewSpXVx = ViewSpXVx<double>;

using DBSViewX = BSViewX<double>;

using RDimFx = Fourier<RDimX>;
using CoordFx = Coordinate<RDimFx>;
using IDimFx = NonUniformPointSampling<RDimFx>;
using IndexFx = DiscreteElement<IDimFx>;
using IVectFx = DiscreteVector<IDimFx>;
using IDomainFx = DiscreteDomain<IDimFx>;
