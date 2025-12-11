

# File geometry.hpp

[**File List**](files.md) **>** [**geometry**](dir_6ef3b5c953c12640e6eb10271de0236d.md) **>** [**geometry.hpp**](geometryXY_2geometry_2geometry_8hpp.md)

[Go to the documentation of this file](geometryXY_2geometry_2geometry_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
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


using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

// IDim definitions
struct GridX : UniformGridBase<X>
{
};
struct GridY : UniformGridBase<Y>
{
};


// Index definitions
using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;

// Index Step definitions
using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;

// Index Range definitions
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
        VectorIndexSet<X, Y>,
        Kokkos::DefaultExecutionSpace::memory_space>;
using VectorFieldXY_XY = typename VectorFieldMemXY_XY::span_type;
using VectorConstFieldXY_XY = typename VectorFieldMemXY_XY::view_type;
```


