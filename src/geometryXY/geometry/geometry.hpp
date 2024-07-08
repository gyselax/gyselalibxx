// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_helper.hpp"
#include "directional_tag.hpp"
#include "vector_field.hpp"
#include "vector_field_span.hpp"



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


using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;
using CoordXY = ddc::Coordinate<RDimX, RDimY>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeY = 3;

bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsY = true;

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

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::PERIODIC;

// IDim initialisers
using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;

// IDim definitions
struct IDimX : SplineInterpPointsX::interpolation_mesh_type
{
};
struct IDimY : SplineInterpPointsY::interpolation_mesh_type
{
};


// SplineBuilder and SplineEvaluator definitions
using SplineXBuilder_XY = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY>;
using SplineXEvaluator_XY = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimX,
        IDimY>;


using SplineYBuilder_XY = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX,
        IDimY>;
using SplineYEvaluator_XY = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        IDimY,
        ddc::PeriodicExtrapolationRule<RDimY>,
        ddc::PeriodicExtrapolationRule<RDimY>,
        IDimX,
        IDimY>;

// Spline domain
using BSDomainX = ddc::DiscreteDomain<BSplinesX>;
using BSDomainY = ddc::DiscreteDomain<BSplinesY>;
using BSDomainXY = ddc::DiscreteDomain<BSplinesX, BSplinesY>;

template <class ElementType>
using BSViewXY = device_t<ddc::ChunkSpan<ElementType const, BSDomainXY>>;
using DBSViewXY = BSViewXY<double>;

// Index definitions
using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;

// IVect definitions
using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;

// IDomain definitions
using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;


// Field definitions
template <class ElementType>
using FieldX = device_t<ddc::Chunk<ElementType, IDomainX>>;
using DFieldX = FieldX<double>;

template <class ElementType>
using FieldY = device_t<ddc::Chunk<ElementType, IDomainY>>;
using DFieldY = FieldY<double>;

template <class ElementType>
using FieldXY = device_t<ddc::Chunk<ElementType, IDomainXY>>;
using DFieldXY = FieldXY<double>;


//  Span definitions
template <class ElementType>
using SpanX = device_t<ddc::ChunkSpan<ElementType, IDomainX>>;
using DSpanX = SpanX<double>;

template <class ElementType>
using SpanY = device_t<ddc::ChunkSpan<ElementType, IDomainY>>;
using DSpanY = SpanY<double>;

template <class ElementType>
using SpanXY = device_t<ddc::ChunkSpan<ElementType, IDomainXY>>;
using DSpanXY = SpanXY<double>;


// View definitions
template <class ElementType>
using ViewX = device_t<ddc::ChunkSpan<ElementType const, IDomainX>>;

template <class ElementType>
using ViewY = device_t<ddc::ChunkSpan<ElementType const, IDomainY>>;

template <class ElementType>
using ViewXY = device_t<ddc::ChunkSpan<ElementType const, IDomainXY>>;
using DViewXY = ViewXY<double>;


// VectorField aliases
using VectorFieldXY_XY = VectorField<
        double,
        IDomainXY,
        NDTag<RDimX, RDimY>,
        ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>;
using VectorSpanXY_XY = typename VectorFieldXY_XY::span_type;
using VectorViewXY_XY = typename VectorFieldXY_XY::view_type;