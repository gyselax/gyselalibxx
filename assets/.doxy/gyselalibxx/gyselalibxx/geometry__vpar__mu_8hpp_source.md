

# File geometry\_vpar\_mu.hpp

[**File List**](files.md) **>** [**geometry**](dir_807bff9d645a62665fbe11aeed095652.md) **>** [**geometry\_vpar\_mu.hpp**](geometry__vpar__mu_8hpp.md)

[Go to the documentation of this file](geometry__vpar__mu_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

/*
 * @file geometry_vpar_mu.hpp
 *
 */

struct Vpar
{
    static bool constexpr PERIODIC = false;
};
struct Mu
{
    static bool constexpr PERIODIC = false;
};

// Coord = position of a coordinate in the vector space
using CoordVpar = Coord<Vpar>;
using CoordMu = Coord<Mu>;


struct GridVpar : UniformGridBase<Vpar>
{
};
struct GridMu : UniformGridBase<Mu>
{
};

// Idx = index of the point in the point sampling
using IdxVpar = Idx<GridVpar>;
using IdxMu = Idx<GridMu>;
using IdxVparMu = Idx<GridVpar, GridMu>;
using IdxSpVparMu = Idx<Species, GridVpar, GridMu>;
using IdxSpVpar = Idx<Species, GridVpar>;

// IdxStep = number of grid points between points in a sampling
using IdxStepVpar = IdxStep<GridVpar>;
using IdxStepMu = IdxStep<GridMu>;
using IdxStepVparMu = IdxStep<GridVpar, GridMu>;
using IdxStepSpVparMu = IdxStep<Species, GridVpar, GridMu>;

// IdxRange = to describe the whole index range (or a sub-index range)
using IdxRangeVpar = IdxRange<GridVpar>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeVparMu = IdxRange<GridVpar, GridMu>;
using IdxRangeSpVparMu = IdxRange<Species, GridVpar, GridMu>;
using IdxRangeSpVpar = IdxRange<Species, GridVpar>;

// template for the fields
template <class ElementType>
using FieldMemVpar = FieldMem<ElementType, IdxRangeVpar>;
using DFieldMemVpar = FieldMemVpar<double>;

template <class ElementType>
using FieldMemMu = FieldMem<ElementType, IdxRangeMu>;
using DFieldMemMu = FieldMemMu<double>;

template <class ElementType>
using FieldMemVparMu = FieldMem<ElementType, IdxRangeVparMu>;
using DFieldMemVparMu = FieldMemVparMu<double>;

template <class ElementType>
using FieldMemSpVparMu = FieldMem<ElementType, IdxRangeSpVparMu>;
using DFieldMemSpVparMu = FieldMemSpVparMu<double>;

template <class ElementType>
using FieldMemSpVpar = FieldMem<ElementType, IdxRangeSpVpar>;
using DFieldMemSpVpar = FieldMemSpVpar<double>;

template <class ElementType>
using FieldVpar = Field<ElementType, IdxRangeVpar>;
using DFieldVpar = FieldVpar<double>;

template <class ElementType>
using FieldMu = Field<ElementType, IdxRangeMu>;
using DFieldMu = FieldMu<double>;

template <class ElementType>
using FieldVparMu = Field<ElementType, IdxRangeVparMu>;
using DFieldVparMu = FieldVparMu<double>;

template <class ElementType>
using FieldSpVparMu = Field<ElementType, IdxRangeSpVparMu>;
using DFieldSpVparMu = FieldSpVparMu<double>;

template <class ElementType>
using FieldSpVpar = Field<ElementType, IdxRangeSpVpar>;
using DFieldSpVpar = FieldSpVpar<double>;

template <class ElementType>
using ConstFieldVpar = ConstField<ElementType, IdxRangeVpar>;
using DConstFieldVpar = ConstFieldVpar<double>;

template <class ElementType>
using ConstFieldMu = ConstField<ElementType, IdxRangeMu>;
using DConstFieldMu = ConstFieldMu<double>;

template <class ElementType>
using ConstFieldVparMu = ConstField<ElementType, IdxRangeVparMu>;
using DConstFieldVparMu = ConstFieldVparMu<double>;

template <class ElementType>
using ConstFieldSpVparMu = ConstField<ElementType, IdxRangeSpVparMu>;
using DConstFieldSpVparMu = ConstFieldSpVparMu<double>;

template <class ElementType>
using ConstFieldSpVpar = ConstField<ElementType, IdxRangeSpVpar>;
using DConstFieldSpVpar = ConstFieldSpVpar<double>;
```


