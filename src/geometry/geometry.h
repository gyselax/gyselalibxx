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

using SplineVxBuilder = SplineBuilder<BSplinesVx, BoundCond::HERMITE, BoundCond::HERMITE>;

using MeshX = typename SplineXBuilder::interpolation_mesh_type;

using MeshVx = typename SplineVxBuilder::interpolation_mesh_type;

using MeshFx = NonUniformMesh<Dim::Fx>;

// Species dimension
class MeshSp
{
public:
    using mcoord_type = MCoord<MeshSp>;

    using rdim_type = void;

    constexpr static std::size_t rank()
    {
        return 1;
    }

    constexpr bool operator==(MeshSp const& other) const
    {
        return true;
    }

    constexpr bool operator!=(MeshSp const& other) const
    {
        return !(*this == other);
    }
};

using MCoordX = MCoord<MeshX>;

using MCoordVx = MCoord<MeshVx>;

using MCoordSp = MCoord<MeshSp>;

using MCoordXVx = MCoord<MeshX, MeshVx>;

using MCoordSpXVx = MCoord<MeshSp, MeshX, MeshVx>;

using MCoordFx = MCoord<MeshFx>;

using MLengthX = MLength<MeshX>;

using MLengthVx = MLength<MeshVx>;

using MLengthXVx = MLength<MeshX, MeshVx>;

using MLengthSpXVx = MLength<MeshSp, MeshX, MeshVx>;

using MLengthFx = MLength<MeshFx>;

using MLengthSp = MLength<MeshSp>;

using UniformMDomainX = ProductMDomain<MeshX>;

using UniformMDomainVx = ProductMDomain<MeshVx>;

using UniformMDomainXVx = ProductMDomain<MeshX, MeshVx>;

using MDomainX = UniformMDomainX;

using MDomainVx = UniformMDomainVx;

using MDomainSp = ProductMDomain<MeshSp>;

using MDomainSpVx = ProductMDomain<MeshSp, MeshVx>;

using MDomainSpXVx = ProductMDomain<MeshSp, MeshX, MeshVx>;

using MDomainFx = ProductMDomain<MeshFx>;

template <class ElementType>
using SpanX = BlockSpan<ElementType, MDomainX>;

using DSpanX = SpanX<double>;

template <class ElementType>
using SpanVx = BlockSpan<ElementType, MDomainVx>;

using DSpanVx = SpanVx<double>;

template <class ElementType>
using SpanSpXVx = BlockSpan<ElementType, MDomainSpXVx>;

using DSpanSpXVx = SpanSpXVx<double>;

template <class ElementType>
using ViewX = BlockSpan<ElementType const, MDomainX>;

using DViewX = ViewX<double>;

template <class ElementType>
using ViewVx = BlockSpan<ElementType const, MDomainVx>;

using DViewVx = ViewVx<double>;

template <class ElementType>
using ViewSp = BlockSpan<ElementType const, MDomainSp>;

using DViewSp = ViewSp<double>;

template <class ElementType>
using ViewSpVx = BlockSpan<ElementType const, MDomainSpVx>;

using DViewSpVx = ViewSpVx<double>;

template <class ElementType>
using ViewSpXVx = BlockSpan<ElementType const, MDomainSpXVx>;

using DViewSpXVx = ViewSpXVx<double>;

template <class ElementType>
using BlockX = Block<ElementType, MDomainX>;

using DBlockX = BlockX<double>;

template <class ElementType>
using BlockVx = Block<ElementType, MDomainVx>;

using DBlockVx = BlockVx<double>;

template <class ElementType>
using BlockSp = Block<ElementType, MDomainSp>;

using DBlockSp = BlockSp<double>;

template <class ElementType>
using BlockSpVx = Block<ElementType, MDomainSpVx>;

using DBlockSpVx = BlockSpVx<double>;

template <class ElementType>
using BlockSpXVx = Block<ElementType, MDomainSpXVx>;

using DBlockSpXVx = BlockSpXVx<double>;
