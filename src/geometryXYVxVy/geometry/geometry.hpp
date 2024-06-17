// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>
#include <ddc/kernels/splines.hpp>

#include <ddc_helper.hpp>
#include <species_info.hpp>

/**
 * @brief A class which describes the real space in the first spatial direction X.
 */
struct RDimX
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = true;
};

/**
 * @brief A class which describes the real space in the second spatial direction Y.
 */
struct RDimY
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = true;
};

/**
 * @brief A class which describes the real space in the second velocity direction X.
 */
struct RDimVx
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

/**
 * @brief A class which describes the real space in the second velocity direction Y.
 */
struct RDimVy
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;
using CoordXY = ddc::Coordinate<RDimX, RDimY>;

using CoordVx = ddc::Coordinate<RDimVx>;
using CoordVy = ddc::Coordinate<RDimVy>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeY = 3;

int constexpr BSDegreeVx = 3;
int constexpr BSDegreeVy = 3;

bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsY = true;

bool constexpr BsplineOnUniformCellsVx = true;
bool constexpr BsplineOnUniformCellsVy = true;

struct BSplinesX
    : std::conditional_t<
              BsplineOnUniformCellsX,
              ddc::UniformBSplines<RDimX, BSDegreeX>,
              ddc::NonUniformBSplines<RDimX, BSDegreeX>>
{
};
struct BSplinesY
    : std::conditional_t<
              BsplineOnUniformCellsY,
              ddc::UniformBSplines<RDimY, BSDegreeY>,
              ddc::NonUniformBSplines<RDimY, BSDegreeY>>
{
};

struct BSplinesVx
    : std::conditional_t<
              BsplineOnUniformCellsVx,
              ddc::UniformBSplines<RDimVx, BSDegreeVx>,
              ddc::NonUniformBSplines<RDimVx, BSDegreeVx>>
{
};
struct BSplinesVy
    : std::conditional_t<
              BsplineOnUniformCellsVy,
              ddc::UniformBSplines<RDimVy, BSDegreeVy>,
              ddc::NonUniformBSplines<RDimVy, BSDegreeVy>>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineVyBoundary = ddc::BoundCond::HERMITE;

// IDim initialisers
using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;
using SplineInterpPointsVy
        = ddc::GrevilleInterpolationPoints<BSplinesVy, SplineVyBoundary, SplineVyBoundary>;

// IDim definition
struct IDimX : SplineInterpPointsX::interpolation_mesh_type
{
};
struct IDimY : SplineInterpPointsY::interpolation_mesh_type
{
};
struct IDimVx : SplineInterpPointsVx::interpolation_mesh_type
{
};
struct IDimVy : SplineInterpPointsVy::interpolation_mesh_type
{
};

// SplineBuilder and SplineEvaluator definition
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineYBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineYEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        ddc::PeriodicExtrapolationRule<RDimY>,
        ddc::PeriodicExtrapolationRule<RDimY>,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        ddc::ConstantExtrapolationRule<RDimVx>,
        ddc::ConstantExtrapolationRule<RDimVx>,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineVyBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVy,
        IDimVy,
        SplineVyBoundary,
        SplineVyBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;
using SplineVyEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVy,
        IDimVy,
        ddc::ConstantExtrapolationRule<RDimVy>,
        ddc::ConstantExtrapolationRule<RDimVy>,
        IDimX,
        IDimY,
        IDimVx,
        IDimVy>;

using SplineVxBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::GINKGO,
        IDimVx>;
using SplineVyBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVy,
        IDimVy,
        SplineVyBoundary,
        SplineVyBoundary,
        ddc::SplineSolver::GINKGO,
        IDimVy>;

using BSDomainX = ddc::DiscreteDomain<BSplinesX>;
using BSDomainY = ddc::DiscreteDomain<BSplinesY>;
using BSDomainXY = ddc::DiscreteDomain<BSplinesX, BSplinesY>;
using BSDomainVx = ddc::DiscreteDomain<BSplinesVx>;
using BSDomainVy = ddc::DiscreteDomain<BSplinesVy>;
using BSDomainVxVy = ddc::DiscreteDomain<BSplinesVx, BSplinesVy>;

template <class ElementType>
using BSViewXY = device_t<ddc::ChunkSpan<ElementType const, BSDomainXY>>;
using DBSViewXY = BSViewXY<double>;

// Index
using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using IndexVx = ddc::DiscreteElement<IDimVx>;
using IndexVy = ddc::DiscreteElement<IDimVy>;
using IndexVxVy = ddc::DiscreteElement<IDimVx, IDimVy>;
using IndexXYVxVy = ddc::DiscreteElement<IDimX, IDimY, IDimVx, IDimVy>;
using IndexSpXYVxVy = ddc::DiscreteElement<IDimSp, IDimX, IDimY, IDimVx, IDimVy>;

// IVect definition
using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;
using IVectVx = ddc::DiscreteVector<IDimVx>;
using IVectVy = ddc::DiscreteVector<IDimVy>;

// Idomain definition
using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainVx = ddc::DiscreteDomain<IDimVx>;
using IDomainVy = ddc::DiscreteDomain<IDimVy>;
using IDomainXYVxVy = ddc::DiscreteDomain<IDimX, IDimY, IDimVx, IDimVy>;
using IDomainVxVy = ddc::DiscreteDomain<IDimVx, IDimVy>;
using IDomainSpVxVy = ddc::DiscreteDomain<IDimSp, IDimVx, IDimVy>;
using IDomainSpXYVxVy = ddc::DiscreteDomain<IDimSp, IDimX, IDimY, IDimVx, IDimVy>;

template <class ElementType>
using FieldX = device_t<ddc::Chunk<ElementType, IDomainX>>;
using DFieldX = FieldX<double>;

template <class ElementType>
using FieldY = device_t<ddc::Chunk<ElementType, IDomainY>>;
using DFieldY = FieldY<double>;

template <class ElementType>
using FieldXY = device_t<ddc::Chunk<ElementType, IDomainXY>>;
using DFieldXY = FieldXY<double>;

template <class ElementType>
using FieldVx = device_t<ddc::Chunk<ElementType, IDomainVx>>;

template <class ElementType>
using FieldVy = device_t<ddc::Chunk<ElementType, IDomainVy>>;

template <class ElementType>
using FieldVxVy = device_t<ddc::Chunk<ElementType, IDomainVxVy>>;
using DFieldVxVy = FieldVxVy<double>;

template <class ElementType>
using FieldXYVxVy = device_t<ddc::Chunk<ElementType, IDomainXYVxVy>>;
using DFieldXYVxVy = FieldXYVxVy<double>;

template <class ElementType>
using FieldSpVxVy = device_t<ddc::Chunk<ElementType, IDomainSpVxVy>>;
using DFieldSpVxVy = FieldSpVxVy<double>;

template <class ElementType>
using FieldSpXYVxVy = device_t<ddc::Chunk<ElementType, IDomainSpXYVxVy>>;
using DFieldSpXYVxVy = FieldSpXYVxVy<double>;

//  Span definition
template <class ElementType>
using SpanX = device_t<ddc::ChunkSpan<ElementType, IDomainX>>;
using DSpanX = SpanX<double>;

template <class ElementType>
using SpanY = device_t<ddc::ChunkSpan<ElementType, IDomainY>>;
using DSpanY = SpanY<double>;

template <class ElementType>
using SpanXY = device_t<ddc::ChunkSpan<ElementType, IDomainXY>>;
using DSpanXY = SpanXY<double>;

template <class ElementType>
using SpanVx = device_t<ddc::ChunkSpan<ElementType, IDomainVx>>;
using DSpanVx = SpanVx<double>;

template <class ElementType>
using SpanVy = device_t<ddc::ChunkSpan<ElementType, IDomainVy>>;
using DSpanVy = SpanVy<double>;

template <class ElementType>
using SpanVxVy = device_t<ddc::ChunkSpan<ElementType, IDomainVxVy>>;
using DSpanVxVy = SpanVxVy<double>;

template <class ElementType>
using SpanSpVxVy = device_t<ddc::ChunkSpan<ElementType, IDomainSpVxVy>>;
using DSpanSpVxVy = SpanSpVxVy<double>;

template <class ElementType>
using SpanSpXYVxVy = device_t<ddc::ChunkSpan<ElementType, IDomainSpXYVxVy>>;
using DSpanSpXYVxVy = SpanSpXYVxVy<double>;

// View definition
template <class ElementType>
using ViewX = device_t<ddc::ChunkSpan<ElementType const, IDomainX>>;

template <class ElementType>
using ViewY = device_t<ddc::ChunkSpan<ElementType const, IDomainY>>;

template <class ElementType>
using ViewXY = device_t<ddc::ChunkSpan<ElementType const, IDomainXY>>;
using DViewXY = ViewXY<double>;

template <class ElementType>
using ViewVx = device_t<ddc::ChunkSpan<ElementType const, IDomainVx>>;

template <class ElementType>
using ViewVy = device_t<ddc::ChunkSpan<ElementType const, IDomainVy>>;

template <class ElementType>
using ViewVxVy = device_t<ddc::ChunkSpan<ElementType const, IDomainVxVy>>;
using DViewVxVy = ViewVxVy<double>;

template <class ElementType>
using ViewSpVxVy = device_t<ddc::ChunkSpan<ElementType const, IDomainSpVxVy>>;
using DViewSpVxVy = ViewSpVxVy<double>;

template <class ElementType>
using ViewSpXYVxVy = device_t<ddc::ChunkSpan<ElementType const, IDomainSpXYVxVy>>;
using DViewSpXYVxVy = ViewSpXYVxVy<double>;

/**
 * @brief A class providing aliases for useful subdomains of the geometry. It is used as template parameter for generic dimensionality-agnostic operat
ors such as advections.
 */
class GeometryXYVxVy
{
public:
    /**
     * @brief A templated type giving the velocity discrete dimension type associated to a spatial discrete dimension type.
     */
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, IDimX>,
            IDimVx,
            std::conditional_t<std::is_same_v<T, IDimY>, IDimVy, void>>;

    /**
     * @brief A templated type giving the spatial discrete dimension type associated to a velocity discrete dimension type.
     */
    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, IDimVx>, IDimX, std::conditional_t<std::is_same_v<T, IDimVy>, IDimY, void>>;

    /**
     * @brief An alias for the spatial discrete domain type.
     */
    using SpatialDDom = IDomainXY;

    /**
     * @brief An alias for the velocity discrete domain type.
     */
    using VelocityDDom = IDomainVxVy;

    /**
     * @brief An alias for the whole distribution function discrete domain type.
     */
    using FdistribuDDom = IDomainSpXYVxVy;
};
