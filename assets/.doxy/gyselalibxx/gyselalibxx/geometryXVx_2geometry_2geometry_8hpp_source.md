

# File geometry.hpp

[**File List**](files.md) **>** [**geometry**](dir_3d8ac113f1c21fd7bc019cd952574dfc.md) **>** [**geometry.hpp**](geometryXVx_2geometry_2geometry_8hpp.md)

[Go to the documentation of this file](geometryXVx_2geometry_2geometry_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "moments.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "species_info.hpp"

struct X
{
#ifdef PERIODIC_RDIMX
    static bool constexpr PERIODIC = true;
#else
    static bool constexpr PERIODIC = false;
#endif
    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = X;
};

struct Vx
{
    static bool constexpr PERIODIC = false;
};


using CoordX = Coord<X>;

using CoordVx = Coord<Vx>;

using CoordXVx = Coord<X, Vx>;

#ifdef INPUT_MESH
bool constexpr GRID_X_UNIFORM = false;
bool constexpr GRID_VX_UNIFORM = false;
#else
bool constexpr GRID_X_UNIFORM = X::PERIODIC;
bool constexpr GRID_VX_UNIFORM = true;
#endif

struct GridX : std::conditional_t<GRID_X_UNIFORM, UniformGridBase<X>, NonUniformGridBase<X>>
{
};
struct GridVx : std::conditional_t<GRID_VX_UNIFORM, UniformGridBase<Vx>, NonUniformGridBase<Vx>>
{
};

struct GridMom : Moments
{
};

using IdxMom = Idx<GridMom>;

using IdxVx = Idx<GridVx>;

using IdxX = Idx<GridX>;


using IdxSpMom = Idx<Species, GridMom>;

using IdxSpMomX = Idx<Species, GridMom, GridX>;

using IdxSpX = Idx<Species, GridX>;

using IdxSpVx = Idx<Species, GridVx>;

using IdxSpXVx = Idx<Species, GridX, GridVx>;

using IdxXVx = Idx<GridX, GridVx>;



using IdxStepMom = IdxStep<GridMom>;

using IdxStepVx = IdxStep<GridVx>;

using IdxStepX = IdxStep<GridX>;


using IdxStepSpMom = IdxStep<Species, GridMom>;

using IdxStepSpMomX = IdxStep<Species, GridMom, GridX>;

using IdxStepSpVx = IdxStep<Species, GridVx>;

using IdxStepSpX = IdxStep<Species, GridX>;

using IdxStepSpXVx = IdxStep<Species, GridX, GridVx>;

using IdxStepXVx = IdxStep<GridX, GridVx>;



using IdxRangeMom = IdxRange<GridMom>;

using IdxRangeVx = IdxRange<GridVx>;

using IdxRangeX = IdxRange<GridX>;

using IdxRangeSpMom = IdxRange<Species, GridMom>;

using IdxRangeSpMomX = IdxRange<Species, GridMom, GridX>;

using IdxRangeSpVx = IdxRange<Species, GridVx>;

using IdxRangeSpX = IdxRange<Species, GridX>;

using IdxRangeSpXVx = IdxRange<Species, GridX, GridVx>;

using IdxRangeXVx = IdxRange<GridX, GridVx>;


template <class ElementType>
using FieldMemVx = FieldMem<ElementType, IdxRangeVx>;

template <class ElementType>
using FieldMemX = FieldMem<ElementType, IdxRangeX>;


template <class ElementType>
using FieldMemSpMom = FieldMem<ElementType, IdxRangeSpMom>;

template <class ElementType>
using FieldMemSpMomX = FieldMem<ElementType, IdxRangeSpMomX>;

template <class ElementType>
using FieldMemSpVx = FieldMem<ElementType, IdxRangeSpVx>;

template <class ElementType>
using FieldMemSpX = FieldMem<ElementType, IdxRangeSpX>;

template <class ElementType>
using FieldMemSpXVx = FieldMem<ElementType, IdxRangeSpXVx>;



using DFieldMemVx = FieldMemVx<double>;

using DFieldMemX = FieldMemX<double>;


using DFieldMemSpMom = FieldMemSpMom<double>;

using DFieldMemSpMomX = FieldMemSpMomX<double>;

using DFieldMemSpVx = FieldMemSpVx<double>;

using DFieldMemSpX = FieldMemSpX<double>;

using DFieldMemSpXVx = FieldMemSpXVx<double>;



template <class ElementType>
using FieldX = Field<ElementType, IdxRangeX>;

template <class ElementType>
using FieldVx = Field<ElementType, IdxRangeVx>;

template <class ElementType>
using FieldSpMomX = Field<ElementType, IdxRangeSpMomX>;

template <class ElementType>
using FieldSpMom = Field<ElementType, IdxRangeSpMom>;

template <class ElementType>
using FieldSpVx = Field<ElementType, IdxRangeSpVx>;

template <class ElementType>
using FieldSpX = Field<ElementType, IdxRangeSpX>;

template <class ElementType>
using FieldSpXVx = Field<ElementType, IdxRangeSpXVx>;


using DFieldVx = FieldVx<double>;

using DFieldX = FieldX<double>;

using DFieldSpMomX = FieldSpMomX<double>;

using DFieldSpMom = FieldSpMom<double>;

using DFieldSpVx = FieldSpVx<double>;

using DFieldSpX = FieldSpX<double>;

using DFieldSpXVx = FieldSpXVx<double>;


template <class ElementType>
using ConstFieldVx = Field<ElementType const, IdxRangeVx>;

template <class ElementType>
using ConstFieldX = Field<ElementType const, IdxRangeX>;


template <class ElementType>
using ConstFieldSpMom = ConstField<ElementType, IdxRangeSpMom>;

template <class ElementType>
using ConstFieldSpMomX = ConstField<ElementType, IdxRangeSpMomX>;

template <class ElementType>
using ConstFieldSpMom = ConstField<ElementType, IdxRangeSpMom>;

template <class ElementType>
using ConstFieldSpVx = ConstField<ElementType, IdxRangeSpVx>;

template <class ElementType>
using ConstFieldSpX = ConstField<ElementType, IdxRangeSpX>;

template <class ElementType>
using ConstFieldSpXVx = ConstField<ElementType, IdxRangeSpXVx>;



using DConstFieldVx = ConstFieldVx<double>;

using DConstFieldX = ConstFieldX<double>;


using DConstFieldSpMom = ConstFieldSpMom<double>;

using DConstFieldSpMomX = ConstFieldSpMomX<double>;

using DConstFieldSpMom = ConstFieldSpMom<double>;

using DConstFieldSpX = ConstFieldSpX<double>;

using DConstFieldSpVx = ConstFieldSpVx<double>;

using DConstFieldSpXVx = ConstFieldSpXVx<double>;


class GeometryXVx
{
public:
    template <class T>
    using velocity_dim_for = std::conditional_t<std::is_same_v<T, GridX>, GridVx, void>;

    template <class T>
    using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, void>;

    using IdxRangeSpatial = IdxRangeX;

    using IdxRangeVelocity = IdxRangeVx;


    // using FdistribuIdxRange = IdxRange<DimSp, typename decltype(SpatialDDom), typename decltype(VelocityDDom)>(IdxRange());
    using IdxRangeFdistribu = IdxRangeSpXVx;
};
```


