// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "geometry_xvx.hpp"
#include "mpilayout.hpp"
#include "species_info.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"

using Real = GYSELALIBXX_BUILD_REAL_PRECISION;

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
 * @brief A class which describes the real space in the second velocity direction Y.
 */
struct Vy
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

using CoordVy = Coord<Vy>;

// IDim definition
struct GridY : UniformGridBase<Y>
{
};
struct GridVy : UniformGridBase<Vy>
{
};

// Index
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;
using IdxVy = Idx<GridVy>;
using IdxVxVy = Idx<GridVx, GridVy>;
using IdxXYVxVy = Idx<GridX, GridY, GridVx, GridVy>;
using IdxSpXYVxVy = Idx<Species, GridX, GridY, GridVx, GridVy>;

// IVect definition
using IdxStepY = IdxStep<GridY>;
using IdxStepVy = IdxStep<GridVy>;

// Iindex range definition
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeVy = IdxRange<GridVy>;
using IdxRangeXYVxVy = IdxRange<GridX, GridY, GridVx, GridVy>;
using IdxRangeVxVyXY = IdxRange<GridVx, GridVy, GridX, GridY>;
using IdxRangeVxVy = IdxRange<GridVx, GridVy>;
using IdxRangeSpVxVy = IdxRange<Species, GridVx, GridVy>;
using IdxRangeSpXYVxVy = IdxRange<Species, GridX, GridY, GridVx, GridVy>;
using IdxRangeSpVxVyXY = IdxRange<Species, GridVx, GridVy, GridX, GridY>;

template <class ElementType>
using FieldMemY = FieldMem<ElementType, IdxRangeY>;
using DFieldMemY = FieldMemY<Real>;

template <class ElementType>
using FieldMemXY = FieldMem<ElementType, IdxRangeXY>;
using DFieldMemXY = FieldMemXY<Real>;

template <class ElementType>
using VectorFieldMemXY = VectorFieldMem<ElementType, IdxRangeXY, VectorIndexSet<X, Y>>;
using DVectorFieldMemXY = VectorFieldMemXY<Real>;

template <class ElementType>
using FieldMemVy = FieldMem<ElementType, IdxRangeVy>;
using DFieldMemVy = FieldMemVy<Real>;

template <class ElementType>
using FieldMemVxVy = FieldMem<ElementType, IdxRangeVxVy>;
using DFieldMemVxVy = FieldMemVxVy<Real>;

template <class ElementType>
using FieldMemXYVxVy = FieldMem<ElementType, IdxRangeXYVxVy>;
using DFieldMemXYVxVy = FieldMemXYVxVy<Real>;

template <class ElementType>
using FieldMemSpVxVy = FieldMem<ElementType, IdxRangeSpVxVy>;
using DFieldMemSpVxVy = FieldMemSpVxVy<Real>;

template <class ElementType>
using FieldMemSpXYVxVy = FieldMem<ElementType, IdxRangeSpXYVxVy>;
using DFieldMemSpXYVxVy = FieldMemSpXYVxVy<Real>;

template <class ElementType>
using FieldMemSpVxVyXY = FieldMem<ElementType, IdxRangeSpVxVyXY>;
using DFieldMemSpVxVyXY = FieldMemSpVxVyXY<Real>;

template <class ElementType>
using FieldY = Field<ElementType, IdxRangeY>;
using DFieldY = FieldY<Real>;

template <class ElementType>
using FieldXY = Field<ElementType, IdxRangeXY>;
using DFieldXY = FieldXY<Real>;

template <class ElementType>
using VectorFieldXY = VectorField<ElementType, IdxRangeXY, VectorIndexSet<X, Y>>;
using DVectorFieldXY = VectorFieldXY<Real>;

template <class ElementType>
using FieldVy = Field<ElementType, IdxRangeVy>;
using DFieldVy = FieldVy<Real>;

template <class ElementType>
using FieldVxVy = Field<ElementType, IdxRangeVxVy>;
using DFieldVxVy = FieldVxVy<Real>;

template <class ElementType>
using FieldSpVxVy = Field<ElementType, IdxRangeSpVxVy>;
using DFieldSpVxVy = FieldSpVxVy<Real>;

template <class ElementType>
using FieldSpXYVxVy = Field<ElementType, IdxRangeSpXYVxVy>;
using DFieldSpXYVxVy = FieldSpXYVxVy<Real>;

template <class ElementType>
using FieldSpVxVyXY = Field<ElementType, IdxRangeSpVxVyXY>;
using DFieldSpVxVyXY = FieldSpVxVyXY<Real>;

template <class ElementType>
using ConstFieldY = Field<ElementType const, IdxRangeY>;

template <class ElementType>
using ConstFieldXY = Field<ElementType const, IdxRangeXY>;
using DConstFieldXY = ConstFieldXY<Real>;

template <class ElementType>
using VectorConstFieldXY = VectorConstField<ElementType, IdxRangeXY, VectorIndexSet<X, Y>>;
using DVectorConstFieldXY = VectorConstFieldXY<Real>;

template <class ElementType>
using ConstFieldVx = Field<ElementType const, IdxRangeVx>;

template <class ElementType>
using ConstFieldVy = Field<ElementType const, IdxRangeVy>;

template <class ElementType>
using ConstFieldVxVy = Field<ElementType const, IdxRangeVxVy>;
using DConstFieldVxVy = ConstFieldVxVy<Real>;

template <class ElementType>
using ConstFieldSpVxVy = Field<ElementType const, IdxRangeSpVxVy>;
using DConstFieldSpVxVy = ConstFieldSpVxVy<Real>;

template <class ElementType>
using ConstFieldSpXYVxVy = Field<ElementType const, IdxRangeSpXYVxVy>;
using DConstFieldSpXYVxVy = ConstFieldSpXYVxVy<Real>;

template <class ElementType>
using ConstFieldSpVxVyXY = Field<ElementType const, IdxRangeSpVxVyXY>;
using DConstFieldSpVxVyXY = ConstFieldSpVxVyXY<Real>;

using X2DSplit = MPILayout<IdxRangeSpXYVxVy, GridX, GridY>;
using V2DSplit = MPILayout<IdxRangeSpVxVyXY, GridVx, GridVy>;

/**
 * @brief A class providing aliases for useful subindex ranges of the geometry when the data is saved with the spatial dimensions
 * distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections.
 */
class GeometryXYVxVy
{
public:
    /**
     * @brief A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type.
     */
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, GridX>,
            GridVx,
            std::conditional_t<std::is_same_v<T, GridY>, GridVy, void>>;

    /**
     * @brief A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type.
     */
    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, std::conditional_t<std::is_same_v<T, GridVy>, GridY, void>>;

    /**
     * @brief An alias for the spatial discrete index range type.
     */
    using IdxRangeSpatial = IdxRangeXY;

    /**
     * @brief An alias for the velocity discrete index range type.
     */
    using IdxRangeVelocity = IdxRangeVxVy;

    /**
     * @brief An alias for the whole distribution function discrete index range type.
     */
    using IdxRangeFdistribu = IdxRangeSpXYVxVy;
};

/**
 * @brief A class providing aliases for useful subindex ranges of the geometry when the data is saved with the velocity dimensions
 * distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections.
 */
class GeometryVxVyXY
{
public:
    /**
     * @brief A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type.
     */
    template <class T>
    using velocity_dim_for = std::conditional_t<
            std::is_same_v<T, GridX>,
            GridVx,
            std::conditional_t<std::is_same_v<T, GridY>, GridVy, void>>;

    /**
     * @brief A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type.
     */
    // template <class T>
    // using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, std::conditional_t<std::is_same_v<T, GridVy>, GridY, void>>;

    /**
     * @brief An alias for the spatial discrete index range type.
     */
    using IdxRangeSpatial = IdxRangeXY;

    /**
     * @brief An alias for the velocity discrete index range type.
     */
    using IdxRangeVelocity = IdxRangeVxVy;

    /**
     * @brief An alias for the whole distribution function discrete index range type.
     */
    using IdxRangeFdistribu = IdxRangeSpVxVyXY;
};
