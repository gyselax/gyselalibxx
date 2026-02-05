// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "mpilayout.hpp"
#include "species_info.hpp"

/**
 * @brief A class which describes the real space in the first spatial direction X.
 */
struct X
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = true;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the co/contra-variant counterpart.
    using Dual = X;
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

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the co/contra-variant counterpart.
    using Dual = Y;
};

/**
 * @brief A class which describes the real space in the second velocity direction X.
 */
struct Vx
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

/**
 * @brief A class which describes the real space in the second velocity direction Y.
 */
struct Vy
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

/**
 * @brief A class which describes the real space in the second velocity direction Y.
 */
struct Vz
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

using CoordX = Coord<X>;

using CoordVx = Coord<Vx>;
using CoordVy = Coord<Vy>;
using CoordVz = Coord<Vz>;

int constexpr BSDegreeX = 3;

int constexpr BSDegreeVx = 3;
int constexpr BSDegreeVy = 3;
int constexpr BSDegreeVz = 3;

bool constexpr BsplineOnUniformCellsX = true;

bool constexpr BsplineOnUniformCellsVx = true;
bool constexpr BsplineOnUniformCellsVy = true;
bool constexpr BsplineOnUniformCellsVz = true;

struct BSplinesX
    : std::conditional_t<
              BsplineOnUniformCellsX,
              ddc::UniformBSplines<X, BSDegreeX>,
              ddc::NonUniformBSplines<X, BSDegreeX>>
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
struct BSplinesVz
    : std::conditional_t<
              BsplineOnUniformCellsVz,
              ddc::UniformBSplines<Vz, BSDegreeVz>,
              ddc::NonUniformBSplines<Vz, BSDegreeVz>>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineVyBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineVzBoundary = ddc::BoundCond::HERMITE;

// IDim initialisers
using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;
using SplineInterpPointsVy
        = ddc::GrevilleInterpolationPoints<BSplinesVy, SplineVyBoundary, SplineVyBoundary>;
using SplineInterpPointsVz
        = ddc::GrevilleInterpolationPoints<BSplinesVz, SplineVyBoundary, SplineVyBoundary>;

// IDim definition
struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridVx : SplineInterpPointsVx::interpolation_discrete_dimension_type
{
};
struct GridVy : SplineInterpPointsVy::interpolation_discrete_dimension_type
{
};
struct GridVz : SplineInterpPointsVz::interpolation_discrete_dimension_type
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
using SplineVzBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVz,
        GridVz,
        SplineVzBoundary,
        SplineVzBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineVzEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVz,
        GridVz,
        ddc::ConstantExtrapolationRule<Vz>,
        ddc::ConstantExtrapolationRule<Vz>>;

using IdxRangeBSX = IdxRange<BSplinesX>;
using IdxRangeBSVx = IdxRange<BSplinesVx>;
using IdxRangeBSVy = IdxRange<BSplinesVy>;
using IdxRangeBSVz = IdxRange<BSplinesVz>;
using IdxRangeBSVxVyVz = IdxRange<BSplinesVx, BSplinesVy, BSplinesVz>;

// Index
using IdxX = Idx<GridX>;
using IdxVx = Idx<GridVx>;
using IdxVy = Idx<GridVy>;
using IdxVz = Idx<GridVz>;

using IdxVxVy = Idx<GridVx, GridVy>;
using IdxVyVz = Idx<GridVy, GridVz>;

using IdxVxVyVz = Idx<GridVx, GridVy, GridVz>;
using IdxXVxVyVz = Idx<GridX, GridVx, GridVy, GridVz>;
using IdxSpXVxVyVz = Idx<Species, GridX, GridVx, GridVy, GridVz>;

// IVect definition
using IdxStepX = IdxStep<GridX>;
using IdxStepVx = IdxStep<GridVx>;
using IdxStepVy = IdxStep<GridVy>;
using IdxStepVz = IdxStep<GridVz>;

// Iindex range definition
using IdxRangeX = IdxRange<GridX>;
using IdxRangeVx = IdxRange<GridVx>;
using IdxRangeVy = IdxRange<GridVy>;
using IdxRangeVz = IdxRange<GridVz>;
using IdxRangeXVxVyVz = IdxRange<GridX, GridVx, GridVy, GridVz>;
using IdxRangeVxVyVzX = IdxRange<GridVx, GridVy, GridVz, GridX>;
using IdxRangeVxVyVz = IdxRange<GridVx, GridVy, GridVz>;
using IdxRangeVxVy = IdxRange<GridVx, GridVy>;
using IdxRangeSpVxVy = IdxRange<Species, GridVx, GridVy>;
using IdxRangeSpXVxVy = IdxRange<Species, GridX, GridVx, GridVy>;
using IdxRangeSpVxVyX = IdxRange<Species, GridVx, GridVy, GridX>;
using IdxRangeSpVxVyVz = IdxRange<Species, GridVx, GridVy, GridVz>;
using IdxRangeSpX = IdxRange<Species, GridX>;
using IdxRangeSpXVxVyVz = IdxRange<Species, GridX, GridVx, GridVy, GridVz>;
using IdxRangeSpVxVyVzX = IdxRange<Species, GridVx, GridVy, GridVz, GridX>;

template <class ElementType>
using FieldMemX = FieldMem<ElementType, IdxRangeX>;
using DFieldMemX = FieldMemX<double>;

template <class ElementType>
using FieldMemVx = FieldMem<ElementType, IdxRangeVx>;
using DFieldMemVx = FieldMemVx<double>;

template <class ElementType>
using FieldMemVy = FieldMem<ElementType, IdxRangeVy>;
using DFieldMemVy = FieldMemVy<double>;

template <class ElementType>
using FieldMemVz = FieldMem<ElementType, IdxRangeVz>;
using DFieldMemVz = FieldMemVz<double>;

template <class ElementType>
using FieldMemVxVyVz = FieldMem<ElementType, IdxRangeVxVyVz>;
using DFieldMemVxVyVz = FieldMemVxVyVz<double>;

template <class ElementType>
using FieldMemXVxVyVz = FieldMem<ElementType, IdxRangeXVxVyVz>;
using DFieldMemXVxVyVz = FieldMemXVxVyVz<double>;

template <class ElementType>
using FieldMemSpVxVyVz = FieldMem<ElementType, IdxRangeSpVxVyVz>;
using DFieldMemSpVxVyVz = FieldMemSpVxVyVz<double>;

template <class ElementType>
using FieldMemSpVxVy = FieldMem<ElementType, IdxRangeSpVxVy>;
using DFieldMemSpVxVy = FieldMemSpVxVy<double>;

template <class ElementType>
using FieldMemSpX = FieldMem<ElementType, IdxRangeSpX>;
using DFieldMemSpX = FieldMemSpX<double>;

template <class ElementType>
using FieldMemSpXVxVyVz = FieldMem<ElementType, IdxRangeSpXVxVyVz>;
using DFieldMemSpXVxVyVz = FieldMemSpXVxVyVz<double>;

template <class ElementType>
using FieldMemSpVxVyVzX = FieldMem<ElementType, IdxRangeSpVxVyVzX>;
using DFieldMemSpVxVyVzX = FieldMemSpVxVyVzX<double>;

template <class ElementType>
using FieldMemSpXVxVy = FieldMem<ElementType, IdxRangeSpXVxVy>;
using DFieldMemSpXVxVy = FieldMemSpXVxVy<double>;

template <class ElementType>
using FieldMemVxVy = FieldMem<ElementType, IdxRangeVxVy>;
using DFieldMemVxVy = FieldMemVxVy<double>;

//  Field definitions
template <class ElementType>
using FieldX = Field<ElementType, IdxRangeX>;
using DFieldX = FieldX<double>;

template <class ElementType>
using FieldVx = Field<ElementType, IdxRangeVx>;
using DFieldVx = FieldVx<double>;

template <class ElementType>
using FieldVy = Field<ElementType, IdxRangeVy>;
using DFieldVy = FieldVy<double>;

template <class ElementType>
using FieldVxVyVz = Field<ElementType, IdxRangeVxVyVz>;
using DFieldVxVyVz = FieldVxVyVz<double>;

template <class ElementType>
using FieldSpVxVyVz = Field<ElementType, IdxRangeSpVxVyVz>;
using DFieldSpVxVyVz = FieldSpVxVyVz<double>;

template <class ElementType>
using FieldSpX = Field<ElementType, IdxRangeSpX>;
using DFieldSpX = FieldSpX<double>;

template <class ElementType>
using FieldSpXVxVyVz = Field<ElementType, IdxRangeSpXVxVyVz>;
using DFieldSpXVxVyVz = FieldSpXVxVyVz<double>;

template <class ElementType>
using FieldSpXVxVy = Field<ElementType, IdxRangeSpXVxVy>;
using DFieldSpXVxVy = FieldSpXVxVy<double>;

template <class ElementType>
using FieldSpVxVyVzX = Field<ElementType, IdxRangeSpVxVyVzX>;
using DFieldSpVxVyVzX = FieldSpVxVyVzX<double>;

template <class ElementType>
using FieldSpVxVyX = Field<ElementType, IdxRangeSpVxVyX>;
using DFieldSpVxVyX = FieldSpVxVyX<double>;

// ConstField definitions
template <class ElementType>
using ConstFieldX = Field<ElementType const, IdxRangeX>;
using DConstFieldX = ConstFieldX<double>;

template <class ElementType>
using ConstFieldVx = Field<ElementType const, IdxRangeVx>;
using DConstFieldVx = ConstFieldVx<double>;

template <class ElementType>
using ConstFieldSpX = Field<ElementType const, IdxRangeSpX>;
using DConstFieldSpX = ConstFieldSpX<double>;

template <class ElementType>
using ConstFieldVy = Field<ElementType const, IdxRangeVy>;
using DConstFieldVy = ConstFieldVy<double>;

template <class ElementType>
using ConstFieldVz = Field<ElementType const, IdxRangeVz>;
using DConstFieldVz = ConstFieldVz<double>;

template <class ElementType>
using ConstFieldVxVyVz = Field<ElementType const, IdxRangeVxVyVz>;
using DConstFieldVxVyVz = ConstFieldVxVyVz<double>;

template <class ElementType>
using ConstFieldSpVxVyVz = Field<ElementType const, IdxRangeSpVxVyVz>;
using DConstFieldSpVxVyVz = ConstFieldSpVxVyVz<double>;

template <class ElementType>
using ConstFieldSpXVxVyVz = Field<ElementType const, IdxRangeSpXVxVyVz>;
using DConstFieldSpXVxVyVz = ConstFieldSpXVxVyVz<double>;

template <class ElementType>
using ConstFieldSpVxVyVzX = Field<ElementType const, IdxRangeSpVxVyVzX>;
using DConstFieldSpVxVyVzX = ConstFieldSpVxVyVzX<double>;

using X1DSplit = MPILayout<IdxRangeSpXVxVyVz, GridX>;
using V3DSplit = MPILayout<IdxRangeSpVxVyVzX, GridVy, GridVz>;

/**
 * @brief A class providing aliases for useful subindex ranges of the geometry when the data is saved with the spatial dimensions
 * distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections.
 */
class GeometryXVxVyVz
{
public:
    /**
     * @brief A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type.
     */
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, GridX>,
            GridVx,
            void>;

    /**
     * @brief A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type.
     */
    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, std::conditional_t<std::is_same_v<T, GridVy>, GridY, void>>;

    /**
     * @brief An alias for the spatial discrete index range type.
     */
    using IdxRangeSpatial = IdxRangeX;
    using IdxRangeSpSpatial = IdxRangeSpX;

    /**
     * @brief An alias for the velocity discrete index range type.
     */
    using IdxRangeVelocity = IdxRangeVxVyVz;

    /**
     * @brief An alias for the whole distribution function discrete index range type.
     */
    using IdxRangeFdistribu = IdxRangeSpXVxVyVz;
};

/**
 * @brief A class providing aliases for useful subindex ranges of the geometry when the data is saved with the velocity dimensions
 * distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections.
 */
class GeometryVxVyVzX
{
public:
    /**
     * @brief A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type.
     */
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, GridX>,
            GridVx,
            void>;

    /**
     * @brief A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type.
     */
    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, std::conditional_t<std::is_same_v<T, GridVy>, GridY, void>>;

    /**
     * @brief An alias for the spatial discrete index range type.
     */
    using IdxRangeSpatial = IdxRangeX;
    using IdxRangeSpSpatial = IdxRangeSpX;

    /**
     * @brief An alias for the velocity discrete index range type.
     */
    using IdxRangeVelocity = IdxRangeVxVyVz;

    /**
     * @brief An alias for the whole distribution function discrete index range type.
     */
    using IdxRangeFdistribu = IdxRangeSpVxVyVzX;
};
