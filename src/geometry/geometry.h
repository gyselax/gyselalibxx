#pragma once

#include <block.h>
#include <bsplines.h>
#include <non_uniform_mesh.h>
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

using MeshX = UniformMesh<Dim::X>;

using MeshVx = UniformMesh<Dim::Vx>;

using MeshXVx = UniformMesh<Dim::X, Dim::Vx>;

using MCoordX = MCoord<Dim::X>;

using MCoordVx = MCoord<Dim::Vx>;

using MCoordXVx = MCoord<Dim::X, Dim::Vx>;

using MCoordFx = MCoord<Dim::Fx>;

using UniformMDomainX = UniformMDomain<Dim::X>;

using UniformMDomainVx = UniformMDomain<Dim::Vx>;

using UniformMDomainXVx = UniformMDomain<Dim::X, Dim::Vx>;

using MDomainX = UniformMDomain<Dim::X>;

using MDomainVx = UniformMDomain<Dim::Vx>;

using MDomainXVx = UniformMDomain<Dim::X, Dim::Vx>;

using MDomainFx = NonUniformMDomain<Dim::Fx>;

using BSplinesX = BSplines<MDomainX, 4>;

using DBlockSpanX = BlockView<MDomainX, double>;

using DBlockSpanVx = BlockView<MDomainVx, double>;

using DBlockSpanXVx = BlockView<MDomainXVx, double>;

using DBlockViewX = BlockView<MDomainX, double const>;

using DBlockViewVx = BlockView<MDomainVx, double const>;

using DBlockViewXVx = BlockView<MDomainXVx, double const>;

template <class ElementType>
using BlockX = Block<UniformMDomain<Dim::X>, ElementType>;

using DBlockX = BlockX<double>;

template <class ElementType>
using BlockVx = Block<UniformMDomain<Dim::Vx>, ElementType>;

using DBlockVx = BlockVx<double>;

template <class ElementType>
using BlockXVx = Block<UniformMDomain<Dim::X, Dim::Vx>, ElementType>;

using DBlockXVx = BlockXVx<double>;
