// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "directional_tag.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



/**
 * @brief A class which describes the real space in the first spatial direction X.
 */
struct X
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = true;
};

/**
 * @brief A class which describes the real space in the second spatial direction Y.
 */
struct Y
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = true;
};


using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeY = 3;

bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsY = true;

struct BSplinesX
    : std::conditional_t<
              BsplineOnUniformCellsX,
              ddc::UniformBSplines<X, BSDegreeX>,
              ddc::NonUniformBSplines<X, BSDegreeX>>
{
};
struct BSplinesY
    : std::conditional_t<
              BsplineOnUniformCellsY,
              ddc::UniformBSplines<Y, BSDegreeY>,
              ddc::NonUniformBSplines<Y, BSDegreeY>>
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
struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};


// SplineBuilder and SplineEvaluator definitions
using SplineXBuilder_XY = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;
using SplineXEvaluator_XY = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>,
        GridX,
        GridY>;


using SplineYBuilder_XY = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;
using SplineYEvaluator_XY = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        ddc::PeriodicExtrapolationRule<Y>,
        ddc::PeriodicExtrapolationRule<Y>,
        GridX,
        GridY>;

// Spline index range
using BSIdxRangeX = IdxRange<BSplinesX>;
using BSIdxRangeY = IdxRange<BSplinesY>;
using BSIdxRangeXY = IdxRange<BSplinesX, BSplinesY>;

template <class ElementType>
using BSConstFieldXY = Field<ElementType const, BSIdxRangeXY>;
using DBSConstFieldXY = BSConstFieldXY<double>;

// Index definitions
using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;

// IVect definitions
using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;

// IDomain definitions
using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;


// Field definitions
template <class ElementType>
using FieldMemX = FieldMem<ElementType, IdxRangeX>;
using DFieldMemX = FieldMemX<double>;

template <class ElementType>
using FieldMemY = FieldMem<ElementType, IdxRangeY>;
using DFieldMemY = FieldMemY<double>;

template <class ElementType>
using FieldMemXY = FieldMem<ElementType, IdxRangeXY>;
using DFieldMemXY = FieldMemXY<double>;


//  Field definitions
template <class ElementType>
using FieldX = Field<ElementType, IdxRangeX>;
using DFieldX = FieldX<double>;

template <class ElementType>
using FieldY = Field<ElementType, IdxRangeY>;
using DFieldY = FieldY<double>;

template <class ElementType>
using FieldXY = Field<ElementType, IdxRangeXY>;
using DFieldXY = FieldXY<double>;


// ConstField definitions
template <class ElementType>
using ConstFieldX = Field<ElementType const, IdxRangeX>;

template <class ElementType>
using ConstFieldY = Field<ElementType const, IdxRangeY>;

template <class ElementType>
using ConstFieldXY = Field<ElementType const, IdxRangeXY>;
using DConstFieldXY = ConstFieldXY<double>;


// VectorFieldMem aliases
// Represent a vector field (v_x, v_y) on indices (x_i, y_j) : (v_x(x_i, y_j), v_y(x_i,y_j))
using VectorFieldMemXY_XY = VectorFieldMem<
        double,
        IdxRangeXY,
        NDTag<X, Y>,
        ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>;
using VectorFieldXY_XY = typename VectorFieldMemXY_XY::span_type;
using VectorConstFieldXY_XY = typename VectorFieldMemXY_XY::view_type;
