#pragma once

#include <ddc/Chunk>
#include <ddc/NonUniformDiscretization>
#include <ddc/UniformDiscretization>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/spline_builder.hpp>

template <class Tag>
struct Fourier
{
    using base_tag_type = Tag;
    static constexpr bool PERIODIC = Tag::PERIODIC;
};



struct RDimX
{
    static constexpr bool PERIODIC = true;
};

struct RDimVx
{
    static constexpr bool PERIODIC = false;
};

struct RDimT
{
    static constexpr bool PERIODIC = false;
};

using RDimFx = Fourier<RDimX>;



using CoordT = Coordinate<RDimT>;

using CoordX = Coordinate<RDimX>;

using CoordVx = Coordinate<RDimVx>;

using CoordXVx = Coordinate<RDimX, RDimVx>;

using CoordFx = Coordinate<RDimFx>;



using BSplinesX = UniformBSplines<RDimX, 3>;

using BSplinesVx = UniformBSplines<RDimVx, 3>;

using SplineXBuilder = SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC>;
using IDimX = typename SplineXBuilder::interpolation_mesh_type;

using SplineVxBuilder = SplineBuilder<BSplinesVx, BoundCond::HERMITE, BoundCond::HERMITE>;
using IDimVx = typename SplineVxBuilder::interpolation_mesh_type;

using IDimFx = NonUniformDiscretization<RDimFx>;

// Species dimension
class IDimSp
{
public:
    using mcoord_type = DiscreteCoordinate<IDimSp>;

    using rdim_type = void;

    constexpr static std::size_t rank()
    {
        return 1;
    }

    constexpr bool operator==(IDimSp const& other) const
    {
        return true;
    }

    constexpr bool operator!=(IDimSp const& other) const
    {
        return !(*this == other);
    }
};



using IndexX = DiscreteCoordinate<IDimX>;

using IndexVx = DiscreteCoordinate<IDimVx>;

using IndexSp = DiscreteCoordinate<IDimSp>;

using IndexXVx = DiscreteCoordinate<IDimX, IDimVx>;

using IndexSpXVx = DiscreteCoordinate<IDimSp, IDimX, IDimVx>;

using IndexFx = DiscreteCoordinate<IDimFx>;



using IVectX = DiscreteVector<IDimX>;

using IVectVx = DiscreteVector<IDimVx>;

using IVectXVx = DiscreteVector<IDimX, IDimVx>;

using IVectSpXVx = DiscreteVector<IDimSp, IDimX, IDimVx>;

using IVectFx = DiscreteVector<IDimFx>;

using IVectSp = DiscreteVector<IDimSp>;



using BSDomainX = DiscreteDomain<BSplinesX>;

using BSDomainVx = DiscreteDomain<BSplinesVx>;

using IDomainX = DiscreteDomain<IDimX>;

using IDomainVx = DiscreteDomain<IDimVx>;

using IDomainXVx = DiscreteDomain<IDimX, IDimVx>;

using IDomainSp = DiscreteDomain<IDimSp>;

using IDomainSpVx = DiscreteDomain<IDimSp, IDimVx>;

using IDomainSpXVx = DiscreteDomain<IDimSp, IDimX, IDimVx>;

using IDomainFx = DiscreteDomain<IDimFx>;



template <class ElementType>
using FieldX = Chunk<ElementType, IDomainX>;

template <class ElementType>
using FieldVx = Chunk<ElementType, IDomainVx>;

template <class ElementType>
using FieldSp = Chunk<ElementType, IDomainSp>;

template <class ElementType>
using FieldSpVx = Chunk<ElementType, IDomainSpVx>;

template <class ElementType>
using FieldSpXVx = Chunk<ElementType, IDomainSpXVx>;



using DFieldX = FieldX<double>;

using DFieldVx = FieldVx<double>;

using DFieldSp = FieldSp<double>;

using DFieldSpVx = FieldSpVx<double>;

using DFieldSpXVx = FieldSpXVx<double>;



template <class ElementType>
using SpanX = ChunkSpan<ElementType, IDomainX>;

template <class ElementType>
using SpanVx = ChunkSpan<ElementType, IDomainVx>;

template <class ElementType>
using SpanSpXVx = ChunkSpan<ElementType, IDomainSpXVx>;

template <class ElementType>
using SpanSpVx = ChunkSpan<ElementType, IDomainSpVx>;



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



using DViewX = ViewX<double>;

using DViewVx = ViewVx<double>;

using DViewSp = ViewSp<double>;

using DViewSpVx = ViewSpVx<double>;

using DViewSpXVx = ViewSpXVx<double>;
