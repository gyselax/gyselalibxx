

# File geometry.hpp

[**File List**](files.md) **>** [**geometry**](dir_7ddd2963f3e4609fce61e92aa9c5ff14.md) **>** [**geometry.hpp**](geometryXYVxVy_2geometry_2geometry_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2geometry_2geometry_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "mpilayout.hpp"
#include "species_info.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"


struct X
{
    static bool constexpr PERIODIC = true;

    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = X;
};

struct Y
{
    static bool constexpr PERIODIC = true;

    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = Y;
};

struct Vx
{
    static bool constexpr PERIODIC = false;
};

struct Vy
{
    static bool constexpr PERIODIC = false;
};

using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

using CoordVx = Coord<Vx>;
using CoordVy = Coord<Vy>;

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

struct BSplinesVx
    : std::conditional_t<
              BsplineOnUniformCellsVx,
              ddc::UniformBSplines<Vx, BSDegreeVx>,
              ddc::NonUniformBSplines<Vx, BSDegreeVx>>
{
};
struct BSplinesVy
    : std::conditional_t<
              BsplineOnUniformCellsVy,
              ddc::UniformBSplines<Vy, BSDegreeVy>,
              ddc::NonUniformBSplines<Vy, BSDegreeVy>>
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
struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};
struct GridVx : SplineInterpPointsVx::interpolation_discrete_dimension_type
{
};
struct GridVy : SplineInterpPointsVy::interpolation_discrete_dimension_type
{
};

// SplineBuilder and SplineEvaluator definition
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>>;
using SplineYBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineYEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        ddc::PeriodicExtrapolationRule<Y>,
        ddc::PeriodicExtrapolationRule<Y>>;
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        ddc::ConstantExtrapolationRule<Vx>,
        ddc::ConstantExtrapolationRule<Vx>>;
using SplineVyBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVy,
        GridVy,
        SplineVyBoundary,
        SplineVyBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineVyEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVy,
        GridVy,
        ddc::ConstantExtrapolationRule<Vy>,
        ddc::ConstantExtrapolationRule<Vy>>;

using IdxRangeBSX = IdxRange<BSplinesX>;
using IdxRangeBSY = IdxRange<BSplinesY>;
using IdxRangeBSXY = IdxRange<BSplinesX, BSplinesY>;
using IdxRangeBSVx = IdxRange<BSplinesVx>;
using IdxRangeBSVy = IdxRange<BSplinesVy>;
using IdxRangeBSVxVy = IdxRange<BSplinesVx, BSplinesVy>;

template <class ElementType>
using BSConstFieldXY = Field<ElementType const, IdxRangeBSXY>;
using DBSConstFieldXY = BSConstFieldXY<double>;

// Index
using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;
using IdxVx = Idx<GridVx>;
using IdxVy = Idx<GridVy>;
using IdxVxVy = Idx<GridVx, GridVy>;
using IdxXYVxVy = Idx<GridX, GridY, GridVx, GridVy>;
using IdxSpXYVxVy = Idx<Species, GridX, GridY, GridVx, GridVy>;

// IVect definition
using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;
using IdxStepVx = IdxStep<GridVx>;
using IdxStepVy = IdxStep<GridVy>;

// Iindex range definition
using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeVx = IdxRange<GridVx>;
using IdxRangeVy = IdxRange<GridVy>;
using IdxRangeXYVxVy = IdxRange<GridX, GridY, GridVx, GridVy>;
using IdxRangeVxVyXY = IdxRange<GridVx, GridVy, GridX, GridY>;
using IdxRangeVxVy = IdxRange<GridVx, GridVy>;
using IdxRangeSpVxVy = IdxRange<Species, GridVx, GridVy>;
using IdxRangeSpXYVxVy = IdxRange<Species, GridX, GridY, GridVx, GridVy>;
using IdxRangeSpVxVyXY = IdxRange<Species, GridVx, GridVy, GridX, GridY>;

template <class ElementType>
using FieldMemX = FieldMem<ElementType, IdxRangeX>;
using DFieldMemX = FieldMemX<double>;

template <class ElementType>
using FieldMemY = FieldMem<ElementType, IdxRangeY>;
using DFieldMemY = FieldMemY<double>;

template <class ElementType>
using FieldMemXY = FieldMem<ElementType, IdxRangeXY>;
using DFieldMemXY = FieldMemXY<double>;

template <class ElementType>
using VectorFieldMemXY = VectorFieldMem<ElementType, IdxRangeXY, VectorIndexSet<X, Y>>;
using DVectorFieldMemXY = VectorFieldMemXY<double>;

template <class ElementType>
using FieldMemVx = FieldMem<ElementType, IdxRangeVx>;

template <class ElementType>
using FieldMemVy = FieldMem<ElementType, IdxRangeVy>;

template <class ElementType>
using FieldMemVxVy = FieldMem<ElementType, IdxRangeVxVy>;
using DFieldMemVxVy = FieldMemVxVy<double>;

template <class ElementType>
using FieldMemXYVxVy = FieldMem<ElementType, IdxRangeXYVxVy>;
using DFieldMemXYVxVy = FieldMemXYVxVy<double>;

template <class ElementType>
using FieldMemSpVxVy = FieldMem<ElementType, IdxRangeSpVxVy>;
using DFieldMemSpVxVy = FieldMemSpVxVy<double>;

template <class ElementType>
using FieldMemSpXYVxVy = FieldMem<ElementType, IdxRangeSpXYVxVy>;
using DFieldMemSpXYVxVy = FieldMemSpXYVxVy<double>;

template <class ElementType>
using FieldMemSpVxVyXY = FieldMem<ElementType, IdxRangeSpVxVyXY>;
using DFieldMemSpVxVyXY = FieldMemSpVxVyXY<double>;

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

template <class ElementType>
using VectorFieldXY = VectorField<ElementType, IdxRangeXY, VectorIndexSet<X, Y>>;
using DVectorFieldXY = VectorFieldXY<double>;

template <class ElementType>
using FieldVx = Field<ElementType, IdxRangeVx>;
using DFieldVx = FieldVx<double>;

template <class ElementType>
using FieldVy = Field<ElementType, IdxRangeVy>;
using DFieldVy = FieldVy<double>;

template <class ElementType>
using FieldVxVy = Field<ElementType, IdxRangeVxVy>;
using DFieldVxVy = FieldVxVy<double>;

template <class ElementType>
using FieldSpVxVy = Field<ElementType, IdxRangeSpVxVy>;
using DFieldSpVxVy = FieldSpVxVy<double>;

template <class ElementType>
using FieldSpXYVxVy = Field<ElementType, IdxRangeSpXYVxVy>;
using DFieldSpXYVxVy = FieldSpXYVxVy<double>;

template <class ElementType>
using FieldSpVxVyXY = Field<ElementType, IdxRangeSpVxVyXY>;
using DFieldSpVxVyXY = FieldSpVxVyXY<double>;

// ConstField definitions
template <class ElementType>
using ConstFieldX = Field<ElementType const, IdxRangeX>;

template <class ElementType>
using ConstFieldY = Field<ElementType const, IdxRangeY>;

template <class ElementType>
using ConstFieldXY = Field<ElementType const, IdxRangeXY>;
using DConstFieldXY = ConstFieldXY<double>;

template <class ElementType>
using VectorConstFieldXY = VectorConstField<ElementType, IdxRangeXY, VectorIndexSet<X, Y>>;
using DVectorConstFieldXY = VectorConstFieldXY<double>;

template <class ElementType>
using ConstFieldVx = Field<ElementType const, IdxRangeVx>;

template <class ElementType>
using ConstFieldVy = Field<ElementType const, IdxRangeVy>;

template <class ElementType>
using ConstFieldVxVy = Field<ElementType const, IdxRangeVxVy>;
using DConstFieldVxVy = ConstFieldVxVy<double>;

template <class ElementType>
using ConstFieldSpVxVy = Field<ElementType const, IdxRangeSpVxVy>;
using DConstFieldSpVxVy = ConstFieldSpVxVy<double>;

template <class ElementType>
using ConstFieldSpXYVxVy = Field<ElementType const, IdxRangeSpXYVxVy>;
using DConstFieldSpXYVxVy = ConstFieldSpXYVxVy<double>;

template <class ElementType>
using ConstFieldSpVxVyXY = Field<ElementType const, IdxRangeSpVxVyXY>;
using DConstFieldSpVxVyXY = ConstFieldSpVxVyXY<double>;

using X2DSplit = MPILayout<IdxRangeSpXYVxVy, GridX, GridY>;
using V2DSplit = MPILayout<IdxRangeSpVxVyXY, GridVx, GridVy>;

class GeometryXYVxVy
{
public:
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, GridX>,
            GridVx,
            std::conditional_t<std::is_same_v<T, GridY>, GridVy, void>>;

    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, std::conditional_t<std::is_same_v<T, GridVy>, GridY, void>>;

    using IdxRangeSpatial = IdxRangeXY;

    using IdxRangeVelocity = IdxRangeVxVy;

    using IdxRangeFdistribu = IdxRangeSpXYVxVy;
};

class GeometryVxVyXY
{
public:
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, GridX>,
            GridVx,
            std::conditional_t<std::is_same_v<T, GridY>, GridVy, void>>;

    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, std::conditional_t<std::is_same_v<T, GridVy>, GridY, void>>;

    using IdxRangeSpatial = IdxRangeXY;

    using IdxRangeVelocity = IdxRangeVxVy;

    using IdxRangeFdistribu = IdxRangeSpVxVyXY;
};
```


