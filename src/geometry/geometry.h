#pragma once

#include <ddc/Block>
#include <ddc/NonUniformMesh>
#include <ddc/UniformMesh>

#include <sll/bsplines.h>
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

using RLengthT = RLength<Dim::T>;

using RLengthX = RLength<Dim::X>;

using RLengthVx = RLength<Dim::Vx>;

using RLengthXVx = RLength<Dim::X, Dim::Vx>;

using KnotsX = UniformMesh<Dim::X>;

using BSplinesX = BSplines<KnotsX, 4>;

using SplineXBuilder = SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC>;

using MeshX = typename SplineXBuilder::interpolation_mesh_type;

using MeshVx = UniformMesh<Dim::Vx>;

using MeshFx = NonUniformMesh<Dim::Fx>;

using MeshXVx = ProductMesh<MeshX, MeshVx>;

using MCoordX = MCoord<MeshX>;

using MCoordVx = MCoord<MeshVx>;

using MCoordXVx = MCoord<MeshX, MeshVx>;

using MCoordFx = MCoord<MeshFx>;

using UniformMDomainX = ProductMDomain<MeshX>;

using UniformMDomainVx = ProductMDomain<MeshVx>;

using UniformMDomainXVx = ProductMDomain<MeshX, MeshVx>;

using MDomainX = UniformMDomainX;

using MDomainVx = UniformMDomainVx;

using MDomainXVx = UniformMDomainXVx;

using MDomainFx = ProductMDomain<MeshFx>;

template <class ElementType>
using SpanX = BlockSpan<MDomainX, ElementType>;

using DSpanX = SpanX<double>;

template <class ElementType>
using SpanVx = BlockSpan<MDomainVx, ElementType>;

using DSpanVx = SpanVx<double>;

template <class ElementType>
using SpanXVx = BlockSpan<MDomainXVx, ElementType>;

using DSpanXVx = SpanXVx<double>;

template <class ElementType>
using ViewX = BlockSpan<MDomainX, ElementType const>;

using DViewX = ViewX<double>;

template <class ElementType>
using ViewVx = BlockSpan<MDomainVx, ElementType const>;

using DViewVx = ViewVx<double>;

template <class ElementType>
using ViewXVx = BlockSpan<MDomainXVx, ElementType const>;

using DViewXVx = ViewXVx<double>;

template <class ElementType>
using BlockX = Block<MDomainX, ElementType>;

using DBlockX = BlockX<double>;

template <class ElementType>
using BlockVx = Block<MDomainVx, ElementType>;

using DBlockVx = BlockVx<double>;

template <class ElementType>
using BlockXVx = Block<MDomainXVx, ElementType>;

using DBlockXVx = BlockXVx<double>;
