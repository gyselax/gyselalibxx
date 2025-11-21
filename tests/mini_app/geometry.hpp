// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "species_info.hpp"

/**
 * @brief Minimal geometry definitions for 5D distribution function
 * This is a minimal version for the mini_app that only defines the necessary types
 * to initialize and read a 5D distribution function (tor1, tor2, tor3, vpar, mu)
 * Note: species is an index, not a dimension
 */

// Toroidal coordinate dimensions
struct Tor1
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Tor1;
};
using Rho = Tor1;

struct Tor2
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Tor2;
};
using Theta = Tor2;

// Velocity space dimensions
struct Vpar
{
    static bool constexpr PERIODIC = false;
};

struct Mu
{
    static bool constexpr PERIODIC = false;
};

// Third toroidal dimension
struct Tor3
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Tor3;
};
using Phi = Tor3;

// Splines definition
int constexpr BSDegreeTor1 = 3;
int constexpr BSDegreeTor2 = 3;
int constexpr BSDegreeTor3 = 3;
int constexpr BSDegreeVpar = 3;
int constexpr BSDegreeMu = 3;

struct BSplinesTor1 : ddc::NonUniformBSplines<Tor1, BSDegreeTor1>
{
};
struct BSplinesTor2 : ddc::NonUniformBSplines<Tor2, BSDegreeTor2>
{
};
struct BSplinesVpar : ddc::NonUniformBSplines<Vpar, BSDegreeVpar>
{
};
struct BSplinesTor3 : ddc::NonUniformBSplines<Tor3, BSDegreeTor3>
{
};
struct BSplinesMu : ddc::NonUniformBSplines<Mu, BSDegreeMu>
{
};

ddc::BoundCond constexpr SplineTor1Boundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineTor2Boundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineTor3Boundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVparBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineMuBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsR = ddcHelper::
        NonUniformInterpolationPoints<BSplinesTor1, SplineTor1Boundary, SplineTor1Boundary>;
using SplineInterpPointsTor2 = ddcHelper::
        NonUniformInterpolationPoints<BSplinesTor2, SplineTor2Boundary, SplineTor2Boundary>;
using SplineInterpPointsTor3 = ddcHelper::
        NonUniformInterpolationPoints<BSplinesTor3, SplineTor3Boundary, SplineTor3Boundary>;
using SplineInterpPointsVpar = ddcHelper::
        NonUniformInterpolationPoints<BSplinesVpar, SplineVparBoundary, SplineVparBoundary>;
using SplineInterpPointsMu
        = ddcHelper::NonUniformInterpolationPoints<BSplinesMu, SplineMuBoundary, SplineMuBoundary>;

struct GridTor1 : NonUniformGridBase<Tor1>
{
};
struct GridTor2 : NonUniformGridBase<Tor2>
{
};
struct GridTor3 : NonUniformGridBase<Tor3>
{
};
struct GridVpar : UniformGridBase<Vpar>
{
};
struct GridMu : UniformGridBase<Mu>
{
};

// IdxRange definitions
using IdxRangeTor1 = IdxRange<GridTor1>;
using IdxRangeTor2 = IdxRange<GridTor2>;
using IdxRangeTor3 = IdxRange<GridTor3>;
using IdxRangeVpar = IdxRange<GridVpar>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeTor2D = IdxRange<GridTor1, GridTor2>;
using IdxRangeTor3D = IdxRange<GridTor1, GridTor2, GridTor3>;
using IdxRangeV2D = IdxRange<GridVpar, GridMu>;
using IdxRangeSpTor3DV2D = IdxRange<Species, GridTor1, GridTor2, GridTor3, GridVpar, GridMu>;

// Field type aliases
template <class ElementType>
using FieldMemSpTor3DV2D = FieldMem<ElementType, IdxRangeSpTor3DV2D>;
using DFieldMemSpTor3DV2D = FieldMemSpTor3DV2D<double>;

template <class ElementType>
using FieldSpTor3DV2D = Field<ElementType, IdxRangeSpTor3DV2D>;
using DFieldSpTor3DV2D = FieldSpTor3DV2D<double>;

template <class ElementType>
using FieldMemSpTor3DV2D_host = host_t<FieldMem<ElementType, IdxRangeSpTor3DV2D>>;
using DFieldMemSpTor3DV2D_host = FieldMemSpTor3DV2D_host<double>;

