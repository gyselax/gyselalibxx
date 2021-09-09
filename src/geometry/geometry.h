#pragma once

#include <block.h>
#include <bsplines.h>
#include <bsplines_non_uniform.h>
#include <bsplines_uniform.h>
#include <non_uniform_mesh.h>
#include <spline_builder.h>
#include <uniform_mesh.h>

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

// using RDomainX = RDomain<Dim::X>;
//
// using RDomainVx = RDomain<Dim::Vx>;
//
// using RDomainXVx = RDomain<Dim::X, Dim::Vx>;

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

using DBlockSpanX = BlockView<MDomainX, double>;

using DBlockSpanVx = BlockView<MDomainVx, double>;

using DBlockSpanXVx = BlockView<MDomainXVx, double>;

using DBlockViewX = BlockView<MDomainX, double const>;

using DBlockViewVx = BlockView<MDomainVx, double const>;

using DBlockViewXVx = BlockView<MDomainXVx, double const>;

template <class ElementType>
using BlockX = Block<UniformMDomainX, ElementType>;

using DBlockX = BlockX<double>;

template <class ElementType>
using BlockVx = Block<UniformMDomainVx, ElementType>;

using DBlockVx = BlockVx<double>;

template <class ElementType>
using BlockXVx = Block<UniformMDomainXVx, ElementType>;

using DBlockXVx = BlockXVx<double>;
