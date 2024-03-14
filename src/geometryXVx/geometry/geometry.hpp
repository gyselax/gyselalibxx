// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>
#include <ddc/kernels/splines.hpp>

#include <ddc_helper.hpp>
#include <species_info.hpp>

/**
 * @brief A class which describes the real space in the spatial X direction
 */
struct RDimX
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
#ifdef PERIODIC_RDIMX
    static bool constexpr PERIODIC = true;
#else
    static bool constexpr PERIODIC = false;
#endif
};

/**
 * @brief A class which describes the real space in the X-velocity direction
 */
struct RDimVx
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

/**
 * @brief A class which describes the real space in the temporal direction
 */
struct RDimT
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};



using CoordT = ddc::Coordinate<RDimT>;

using CoordX = ddc::Coordinate<RDimX>;

using CoordVx = ddc::Coordinate<RDimVx>;

using CoordXVx = ddc::Coordinate<RDimX, RDimVx>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsVx = true;

using BSplinesX = std::conditional_t<
        BsplineOnUniformCellsX,
        ddc::UniformBSplines<RDimX, BSDegreeX>,
        ddc::NonUniformBSplines<RDimX, BSDegreeX>>;
using BSplinesVx = std::conditional_t<
        BsplineOnUniformCellsVx,
        ddc::UniformBSplines<RDimVx, BSDegreeVx>,
        ddc::NonUniformBSplines<RDimVx, BSDegreeVx>>;

auto constexpr SplineXBoundary
        = RDimX::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;
auto constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;

bool constexpr UniformMeshX = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsX,
        SplineXBoundary,
        SplineXBoundary,
        BSDegreeX);
bool constexpr UniformMeshVx = ddc::is_spline_interpolation_mesh_uniform(
        BsplineOnUniformCellsVx,
        SplineVxBoundary,
        SplineVxBoundary,
        BSDegreeVx);

using IDimX = std::conditional_t<
        UniformMeshX,
        ddc::UniformPointSampling<RDimX>,
        ddc::NonUniformPointSampling<RDimX>>;
using IDimVx = std::conditional_t<
        UniformMeshVx,
        ddc::UniformPointSampling<RDimVx>,
        ddc::NonUniformPointSampling<RDimVx>>;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;

using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimVx>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
#ifdef PERIODIC_RDIMX
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
#else
        ddc::ConstantExtrapolationRule<RDimX>,
        ddc::ConstantExtrapolationRule<RDimX>,
#endif
        IDimX,
        IDimVx>;
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimVx>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        ddc::ConstantExtrapolationRule<RDimVx>,
        ddc::ConstantExtrapolationRule<RDimVx>,
        IDimX,
        IDimVx>;
using SplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX>;
using SplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
#ifdef PERIODIC_RDIMX
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
#else
        ddc::ConstantExtrapolationRule<RDimX>,
        ddc::ConstantExtrapolationRule<RDimX>,
#endif
        IDimX>;
using SplineVxBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::GINKGO,
        IDimVx>;
using SplineVxEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        ddc::ConstantExtrapolationRule<RDimVx>,
        ddc::ConstantExtrapolationRule<RDimVx>,
        IDimVx>;

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
using FieldX = device_t<ddc::Chunk<ElementType, IDomainX>>;

template <class ElementType>
using FieldVx = device_t<ddc::Chunk<ElementType, IDomainVx>>;

template <class ElementType>
using FieldSp = device_t<ddc::Chunk<ElementType, IDomainSp>>;

template <class ElementType>
using FieldSpX = device_t<ddc::Chunk<ElementType, IDomainSpX>>;

template <class ElementType>
using FieldSpVx = device_t<ddc::Chunk<ElementType, IDomainSpVx>>;

template <class ElementType>
using FieldSpXVx = device_t<ddc::Chunk<ElementType, IDomainSpXVx>>;



template <class DomainType>
using DField = device_t<ddc::Chunk<double, DomainType>>;


using DFieldX = FieldX<double>;

using DFieldVx = FieldVx<double>;

using DFieldSp = FieldSp<double>;

using DFieldSpX = FieldSpX<double>;

using DFieldSpVx = FieldSpVx<double>;

using DFieldSpXVx = FieldSpXVx<double>;


template <class ElementType>
using SpanSp = device_t<ddc::ChunkSpan<ElementType, IDomainSp>>;

template <class ElementType>
using SpanSpX = device_t<ddc::ChunkSpan<ElementType, IDomainSpX>>;

template <class ElementType>
using SpanX = device_t<ddc::ChunkSpan<ElementType, IDomainX>>;

template <class ElementType>
using SpanSpX = device_t<ddc::ChunkSpan<ElementType, IDomainSpX>>;

template <class ElementType>
using SpanVx = device_t<ddc::ChunkSpan<ElementType, IDomainVx>>;

template <class ElementType>
using SpanSpXVx = device_t<ddc::ChunkSpan<ElementType, IDomainSpXVx>>;

template <class ElementType>
using SpanSpVx = device_t<ddc::ChunkSpan<ElementType, IDomainSpVx>>;

template <class ElementType>
using BSSpanX = device_t<ddc::ChunkSpan<ElementType, BSDomainX>>;



template <class DomainType>
using DSpan = device_t<ddc::ChunkSpan<double, DomainType>>;


using DSpanSp = SpanSp<double>;

using DSpanSpX = SpanSpX<double>;

using DSpanX = SpanX<double>;

using DSpanVx = SpanVx<double>;

using DSpanSpXVx = SpanSpXVx<double>;

using DSpanSpVx = SpanSpVx<double>;

using DBSSpanX = BSSpanX<double>;


template <class ElementType>
using ViewX = device_t<ddc::ChunkSpan<ElementType const, IDomainX>>;

template <class ElementType>
using ViewVx = device_t<ddc::ChunkSpan<ElementType const, IDomainVx>>;

template <class ElementType>
using ViewSp = device_t<ddc::ChunkSpan<ElementType const, IDomainSp>>;

template <class ElementType>
using ViewSpX = device_t<ddc::ChunkSpan<ElementType const, IDomainSpX>>;

template <class ElementType>
using ViewSpVx = device_t<ddc::ChunkSpan<ElementType const, IDomainSpVx>>;

template <class ElementType>
using ViewSpXVx = device_t<ddc::ChunkSpan<ElementType const, IDomainSpXVx>>;

template <class ElementType>
using BSViewX = device_t<ddc::ChunkSpan<ElementType const, BSDomainX>>;



template <class DomainType>
using DView = device_t<ddc::ChunkSpan<double const, DomainType>>;


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

/**
 * @brief A class providing aliases for useful subdomains of the geometry. It is used as template parameter for generic dimensionality-agnostic operators such as advections.
 */
class GeometryXVx
{
public:
    /**
     * @brief A templated type giving the velocity discrete dimension type associated to a spatial discrete dimension type.
     */
    template <class T>
    using velocity_dim_for = std::conditional_t<std::is_same_v<T, IDimX>, IDimVx, void>;

    /**
     * @brief A templated type giving the spatial discrete dimension type associated to a velocity discrete dimension type.
     */
    template <class T>
    using spatial_dim_for = std::conditional_t<std::is_same_v<T, IDimVx>, IDimX, void>;

    /**
     * @brief An alias for species "discrete dimension" type.
     */
    using DDimSp = IDimSp;

    /**
     * @brief An alias for the spatial discrete domain type.
     */
    using SpatialDDom = IDomainX;

    /**
     * @brief An alias for the velocity discrete domain type.
     */
    using VelocityDDom = IDomainVx;


    // using FdistribuDDom = DiscreteDomain<DimSp, typename decltype(SpatialDDom), typename decltype(VelocityDDom)>(ddc::DiscreteDomain());
    /**
     * @brief An alias for the whole distribution function discrete domain type.
     */
    using FdistribuDDom = IDomainSpXVx;
};
