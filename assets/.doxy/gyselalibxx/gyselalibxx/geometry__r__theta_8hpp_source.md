

# File geometry\_r\_theta.hpp

[**File List**](files.md) **>** [**geometry**](dir_718520565cc7a7cfd9ba0e7c9c4c6d52.md) **>** [**geometry\_r\_theta.hpp**](geometry__r__theta_8hpp.md)

[Go to the documentation of this file](geometry__r__theta_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "polar_bsplines.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"


/*
 * @file geometry_r_theta.hpp
 *
 * Definition of
 *   - @f$ r@f$, @f$ \theta@f$, @f$(r, \theta)@f$ dimensions.
 *   - @f$x@f$, @f$y@f$, @f$(x, y)@f$ dimensions.
 */


// POLAR SPACE AND VELOCITY ----------------------------------------------------------------------
// --- Continuous dimensions
struct R_cov;
struct Theta_cov;
struct R
{
    static bool constexpr PERIODIC = false;

    static bool constexpr IS_COVARIANT = false;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = R_cov;
};
struct Theta
{
    static bool constexpr PERIODIC = true;

    static bool constexpr IS_COVARIANT = false;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = Theta_cov;
};

struct R_cov
{
    static bool constexpr PERIODIC = false;

    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = false;

    using Dual = R;
};
struct Theta_cov
{
    static bool constexpr PERIODIC = true;

    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = false;

    using Dual = Theta;
};

struct Vr
{
    static bool constexpr PERIODIC = false;
};
struct Vtheta
{
    static bool constexpr PERIODIC = true;
};


using CoordR = Coord<R>;
using CoordTheta = Coord<Theta>;
using CoordRTheta = Coord<R, Theta>;

using CoordVr = Coord<Vr>;
using CoordVtheta = Coord<Vtheta>;

// --- Discrete dimensions
struct GridR : NonUniformGridBase<R>
{
};
struct GridTheta : NonUniformGridBase<Theta>
{
};
// --- Index definitions
using IdxR = Idx<GridR>;
using IdxTheta = Idx<GridTheta>;
using IdxRTheta = Idx<GridR, GridTheta>;

// --- Index Step definitions
using IdxStepR = IdxStep<GridR>;
using IdxStepTheta = IdxStep<GridTheta>;
using IdxStepRTheta = IdxStep<GridR, GridTheta>;

// --- Index Range definitions
using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;


// --- FieldMem definitions
template <class ElementType>
using FieldMemR = FieldMem<ElementType, IdxRangeR>;

template <class ElementType>
using FieldMemTheta = FieldMem<ElementType, IdxRangeTheta>;

template <class ElementType>
using FieldMemRTheta = FieldMem<ElementType, IdxRangeRTheta>;

using DFieldMemR = FieldMemR<double>;
using DFieldMemTheta = FieldMemTheta<double>;
using DFieldMemRTheta = FieldMemRTheta<double>;

// --- Field definitions
template <class ElementType>
using FieldR = Field<ElementType, IdxRangeR>;

template <class ElementType>
using FieldTheta = Field<ElementType, IdxRangeTheta>;

template <class ElementType>
using FieldRTheta = Field<ElementType, IdxRangeRTheta>;

using DFieldR = FieldR<double>;
using DFieldTheta = FieldTheta<double>;
using DFieldRTheta = FieldRTheta<double>;

// --- Const Field definitions
template <class ElementType>
using ConstFieldR = ConstField<ElementType, IdxRangeR>;

template <class ElementType>
using ConstFieldTheta = ConstField<ElementType, IdxRangeTheta>;

template <class ElementType>
using ConstFieldRTheta = ConstField<ElementType, IdxRangeRTheta>;

using DConstFieldR = ConstFieldR<double>;
using DConstFieldTheta = ConstFieldTheta<double>;
using DConstFieldRTheta = ConstFieldRTheta<double>;


// --- VectorFieldMem definitions
template <class Dim1, class Dim2>
using DVectorFieldMemRTheta = VectorFieldMem<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using DVectorFieldRTheta = VectorField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>>;

template <class Dim1, class Dim2>
using DVectorConstFieldRTheta
        = VectorConstField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>>;



// CARTESIAN SPACE AND VELOCITY ------------------------------------------------------------------
// --- Continuous dimensions
struct X
{
    static bool constexpr PERIODIC = false;

    static bool constexpr IS_COVARIANT = true;

    static bool constexpr IS_CONTRAVARIANT = true;

    using Dual = X;
};
struct Y
{
    static bool constexpr PERIODIC = false;

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
```


