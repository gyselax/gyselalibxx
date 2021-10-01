#pragma once

#include <ddc/Block>
#include <ddc/NonUniformMesh>
#include <ddc/UniformMesh>

#include <sll/bsplines_non_uniform.h>
#include <sll/bsplines_uniform.h>
#include <sll/spline_builder.h>

template <class Tag>
struct Fourier
{
    using base_tag_type = Tag;
    static constexpr bool PERIODIC = Tag::PERIODIC;
};

namespace Dim {

struct X
{
    static constexpr bool PERIODIC = true;
};

struct Vx
{
    static constexpr bool PERIODIC = false;
};

struct T
{
    static constexpr bool PERIODIC = false;
};

using Fx = Fourier<X>;

using Fvx = Fourier<Vx>;

} // namespace Dim

using RCoordT = RCoord<Dim::T>;

using RCoordX = RCoord<Dim::X>;

using RCoordVx = RCoord<Dim::Vx>;

using RCoordXVx = RCoord<Dim::X, Dim::Vx>;

using RCoordFx = RCoord<Dim::Fx>;

using BSplinesX = UniformBSplines<Dim::X, 3>;

using BSDomainX = ProductMDomain<BSplinesX>;

using SplineXBuilder = SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC>;

using BSplinesVx = UniformBSplines<Dim::Vx, 3>;

using BSDomainVx = ProductMDomain<BSplinesVx>;

using SplineVxBuilder = SplineBuilder<BSplinesVx, BoundCond::GREVILLE, BoundCond::GREVILLE>;

using MeshX = typename SplineXBuilder::interpolation_mesh_type;

using MeshVx = UniformMesh<Dim::Vx>;

using MeshFx = NonUniformMesh<Dim::Fx>;

using MeshXVx = ProductMesh<MeshX, MeshVx>;

using MCoordX = MCoord<MeshX>;

using MCoordVx = MCoord<MeshVx>;

using MCoordXVx = MCoord<MeshX, MeshVx>;

using MCoordFx = MCoord<MeshFx>;

using MLengthX = MLength<MeshX>;

using MLengthVx = MLength<MeshVx>;

using MLengthXVx = MLength<MeshX, MeshVx>;

using MLengthFx = MLength<MeshFx>;

using UniformMDomainX = ProductMDomain<MeshX>;

using UniformMDomainVx = ProductMDomain<MeshVx>;

using UniformMDomainXVx = ProductMDomain<MeshX, MeshVx>;

using MDomainX = UniformMDomainX;

using MDomainVx = UniformMDomainVx;

using MDomainXVx = UniformMDomainXVx;

using MDomainFx = ProductMDomain<MeshFx>;

template <class ElementType>
using SpanX = BlockSpan<ElementType, MDomainX>;

using DSpanX = SpanX<double>;

template <class ElementType>
using SpanVx = BlockSpan<ElementType, MDomainVx>;

using DSpanVx = SpanVx<double>;

template <class ElementType>
using SpanXVx = BlockSpan<ElementType, MDomainXVx>;

using DSpanXVx = SpanXVx<double>;

template <class ElementType>
using ViewX = BlockSpan<ElementType const, MDomainX>;

using DViewX = ViewX<double>;

template <class ElementType>
using ViewVx = BlockSpan<ElementType const, MDomainVx>;

using DViewVx = ViewVx<double>;

template <class ElementType>
using ViewXVx = BlockSpan<ElementType const, MDomainXVx>;

using DViewXVx = ViewXVx<double>;

template <class ElementType>
using BlockX = Block<ElementType, MDomainX>;

using DBlockX = BlockX<double>;

template <class ElementType>
using BlockVx = Block<ElementType, MDomainVx>;

using DBlockVx = BlockVx<double>;

template <class ElementType>
using BlockXVx = Block<ElementType, MDomainXVx>;

using DBlockXVx = BlockXVx<double>;
