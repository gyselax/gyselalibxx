// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp> // Strange but needed to get access to Fourier tag, to be clarified

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/spline_boundary_conditions.hpp>
#include <sll/spline_builder.hpp>

#include <species_info.hpp>

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



using CoordT = ddc::Coordinate<RDimT>;

using CoordX = ddc::Coordinate<RDimX>;

using CoordVx = ddc::Coordinate<RDimVx>;

using CoordXVx = ddc::Coordinate<RDimX, RDimVx>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

using BSplinesX = UniformBSplines<RDimX, BSDegreeX>;

using BSplinesVx = UniformBSplines<RDimVx, BSDegreeVx>;

auto constexpr SplineXBoundary = RDimX::PERIODIC ? BoundCond::PERIODIC : BoundCond::GREVILLE;
using InterpPointsX = GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using IDimX = typename InterpPointsX::interpolation_mesh_type;
using SplineXBuilder = SplineBuilder<BSplinesX, IDimX, SplineXBoundary, SplineXBoundary>;

using InterpPointsVx
        = GrevilleInterpolationPoints<BSplinesVx, BoundCond::HERMITE, BoundCond::HERMITE>;
using IDimVx = typename InterpPointsVx::interpolation_mesh_type;
using SplineVxBuilder = SplineBuilder<BSplinesVx, IDimVx, BoundCond::HERMITE, BoundCond::HERMITE>;

// Species dimension
using IDimSp = SpeciesInformation;


using IndexX = ddc::DiscreteElement<IDimX>;

using IndexVx = ddc::DiscreteElement<IDimVx>;

using IndexSp = ddc::DiscreteElement<IDimSp>;

using IndexSpX = ddc::DiscreteElement<IDimSp, IDimX>;

using IndexXVx = ddc::DiscreteElement<IDimX, IDimVx>;

using IndexSpVx = ddc::DiscreteElement<IDimSp, IDimVx>;

using IndexSpXVx = ddc::DiscreteElement<IDimSp, IDimX, IDimVx>;



using IVectX = ddc::DiscreteVector<IDimX>;

using IVectVx = ddc::DiscreteVector<IDimVx>;

using IVectXVx = ddc::DiscreteVector<IDimX, IDimVx>;

using IVectSpXVx = ddc::DiscreteVector<IDimSp, IDimX, IDimVx>;

using IVectSp = ddc::DiscreteVector<IDimSp>;

using IVectSpX = ddc::DiscreteVector<IDimSp, IDimX>;

using IVectSpVx = ddc::DiscreteVector<IDimSp, IDimVx>;



using BSDomainX = ddc::DiscreteDomain<BSplinesX>;

using BSDomainVx = ddc::DiscreteDomain<BSplinesVx>;

using IDomainX = ddc::DiscreteDomain<IDimX>;

using IDomainVx = ddc::DiscreteDomain<IDimVx>;

using IDomainXVx = ddc::DiscreteDomain<IDimX, IDimVx>;

using IDomainSp = ddc::DiscreteDomain<IDimSp>;

using IDomainSpX = ddc::DiscreteDomain<IDimSp, IDimX>;

using IDomainSpVx = ddc::DiscreteDomain<IDimSp, IDimVx>;

using IDomainSpXVx = ddc::DiscreteDomain<IDimSp, IDimX, IDimVx>;



template <class ElementType>
using FieldX = ddc::Chunk<ElementType, IDomainX>;

template <class ElementType>
using FieldVx = ddc::Chunk<ElementType, IDomainVx>;

template <class ElementType>
using FieldSp = ddc::Chunk<ElementType, IDomainSp>;

template <class ElementType>
using FieldSpX = ddc::Chunk<ElementType, IDomainSpX>;

template <class ElementType>
using FieldSpVx = ddc::Chunk<ElementType, IDomainSpVx>;

template <class ElementType>
using FieldSpXVx = ddc::Chunk<ElementType, IDomainSpXVx>;



template <class DomainType>
using DField = ddc::Chunk<double, DomainType>;


using DFieldX = FieldX<double>;

using DFieldVx = FieldVx<double>;

using DFieldSp = FieldSp<double>;

using DFieldSpX = FieldSpX<double>;

using DFieldSpVx = FieldSpVx<double>;

using DFieldSpXVx = FieldSpXVx<double>;


template <class ElementType>
using SpanSp = ddc::ChunkSpan<ElementType, IDomainSp>;

template <class ElementType>
using SpanSpX = ddc::ChunkSpan<ElementType, IDomainSpX>;

template <class ElementType>
using SpanX = ddc::ChunkSpan<ElementType, IDomainX>;

template <class ElementType>
using SpanSpX = ddc::ChunkSpan<ElementType, IDomainSpX>;

template <class ElementType>
using SpanVx = ddc::ChunkSpan<ElementType, IDomainVx>;

template <class ElementType>
using SpanSpXVx = ddc::ChunkSpan<ElementType, IDomainSpXVx>;

template <class ElementType>
using SpanSpVx = ddc::ChunkSpan<ElementType, IDomainSpVx>;



template <class DomainType>
using DSpan = ddc::ChunkSpan<double, DomainType>;


using DSpanSp = SpanSp<double>;

using DSpanSpX = SpanSpX<double>;

using DSpanX = SpanX<double>;

using DSpanVx = SpanVx<double>;

using DSpanSpXVx = SpanSpXVx<double>;

using DSpanSpVx = SpanSpVx<double>;



template <class ElementType>
using ViewX = ddc::ChunkSpan<ElementType const, IDomainX>;

template <class ElementType>
using ViewVx = ddc::ChunkSpan<ElementType const, IDomainVx>;

template <class ElementType>
using ViewSp = ddc::ChunkSpan<ElementType const, IDomainSp>;

template <class ElementType>
using ViewSpX = ddc::ChunkSpan<ElementType const, IDomainSpX>;

template <class ElementType>
using ViewSpVx = ddc::ChunkSpan<ElementType const, IDomainSpVx>;

template <class ElementType>
using ViewSpXVx = ddc::ChunkSpan<ElementType const, IDomainSpXVx>;

template <class ElementType>
using BSViewX = ddc::ChunkSpan<ElementType const, BSDomainX>;



template <class DomainType>
using DView = ddc::ChunkSpan<double const, DomainType>;


using DViewX = ViewX<double>;

using DViewVx = ViewVx<double>;

using DViewSp = ViewSp<double>;

using DViewSpX = ViewSpX<double>;

using DViewSpVx = ViewSpVx<double>;

using DViewSpXVx = ViewSpXVx<double>;

using DBSViewX = BSViewX<double>;

using RDimFx = ddc::Fourier<RDimX>;
using CoordFx = ddc::Coordinate<RDimFx>;
using IDimFx = ddc::PeriodicSampling<RDimFx>;
using IndexFx = ddc::DiscreteElement<IDimFx>;
using IVectFx = ddc::DiscreteVector<IDimFx>;
using IDomainFx = ddc::DiscreteDomain<IDimFx>;

class GeometryXVx
{
public:
    template <class T>
    using velocity_dim_for = std::conditional_t<std::is_same_v<T, IDimX>, IDimVx, void>;

    template <class T>
    using spatial_dim_for = std::conditional_t<std::is_same_v<T, IDimVx>, IDimX, void>;

    using DDimSp = IDimSp;

    using SpatialDDom = IDomainX;

    using VelocityDDom = IDomainVx;


    // using FdistribuDDom = DiscreteDomain<DimSp, typename decltype(SpatialDDom), typename decltype(VelocityDDom)>(ddc::DiscreteDomain());
    using FdistribuDDom = IDomainSpXVx;
};
