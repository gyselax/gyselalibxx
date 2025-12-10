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
 * @file geometry.hpp
 *
 * Definition of
 *   - @f$ r@f$, @f$ \theta@f$, @f$(r, \theta)@f$ dimensions.
 *   - @f$x@f$, @f$y@f$, @f$(x, y)@f$ dimensions.
 */


// POLAR SPACE AND VELOCITY ----------------------------------------------------------------------
// --- Continuous dimensions
struct R_cov;
struct Theta_cov;
/**
 * @brief Define non periodic real contravariant R dimension.
 */
struct R
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = false;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = R_cov;
};
/**
 * @brief Define periodic real contravariant Theta dimension.
 */
struct Theta
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = false;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = Theta_cov;
};

/**
 * @brief Define non periodic real covariant R dimension.
 */
struct R_cov
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = false;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = R;
};
/**
 * @brief Define periodic real covariant Theta dimension.
 */
struct Theta_cov
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = false;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = Theta;
};

/**
 * @brief Define non periodic real R velocity dimension.
 */
struct Vr
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define periodic real Theta velocity dimension.
 */
struct Vtheta
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
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
/**
 * @brief Define non periodic real X dimension.
 */
struct X
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = X;
};
/**
 * @brief Define non periodic real Y dimension.
 */
struct Y
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;

    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;

    /// A type-alias mapping to the covariant counterpart.
    using Dual = Y;
};

/**
 * @brief Define non periodic real X velocity dimension.
 */
struct Vx
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define non periodic real Y velocity dimension.
 */
struct Vy
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};


using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

using CoordVx = Coord<Vx>;
using CoordVy = Coord<Vy>;
