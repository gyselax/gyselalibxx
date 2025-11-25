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
 * @brief Minimal geometry definitions for 6D distribution function
 * This is a minimal version for the mini_app that only defines the necessary types
 * to initialize and read a 6D distribution function (species, tor1, tor2, vpar, mu, dim6)
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

// Additional 6th dimension (not tor3)
struct Dim6
{
    static bool constexpr PERIODIC = false;
};

// Splines definition
int constexpr BSDegreeTor1 = 3;
int constexpr BSDegreeTor2 = 3;
int constexpr BSDegreeVpar = 3;
int constexpr BSDegreeMu = 3;
int constexpr BSDegreeDim6 = 3;

struct BSplinesTor1 : ddc::NonUniformBSplines<Tor1, BSDegreeTor1>
{
};
struct BSplinesTor2 : ddc::NonUniformBSplines<Tor2, BSDegreeTor2>
{
};
struct BSplinesVpar : ddc::NonUniformBSplines<Vpar, BSDegreeVpar>
{
};
struct BSplinesMu : ddc::NonUniformBSplines<Mu, BSDegreeMu>
{
};
struct BSplinesDim6 : ddc::NonUniformBSplines<Dim6, BSDegreeDim6>
{
};

ddc::BoundCond constexpr SplineTor1Boundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineTor2Boundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVparBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineMuBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineDim6Boundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsR = ddcHelper::
        NonUniformInterpolationPoints<BSplinesTor1, SplineTor1Boundary, SplineTor1Boundary>;
using SplineInterpPointsTor2 = ddcHelper::
        NonUniformInterpolationPoints<BSplinesTor2, SplineTor2Boundary, SplineTor2Boundary>;
using SplineInterpPointsVpar = ddcHelper::
        NonUniformInterpolationPoints<BSplinesVpar, SplineVparBoundary, SplineVparBoundary>;
using SplineInterpPointsMu
        = ddcHelper::NonUniformInterpolationPoints<BSplinesMu, SplineMuBoundary, SplineMuBoundary>;
using SplineInterpPointsDim6
        = ddcHelper::NonUniformInterpolationPoints<BSplinesDim6, SplineDim6Boundary, SplineDim6Boundary>;

struct GridTor1 : NonUniformGridBase<Tor1>
{
};
struct GridTor2 : NonUniformGridBase<Tor2>
{
};
struct GridVpar : UniformGridBase<Vpar>
{
};
struct GridMu : UniformGridBase<Mu>
{
};
struct GridDim6 : UniformGridBase<Dim6>
{
};

// IdxRange definitions
using IdxRangeTor1 = IdxRange<GridTor1>;
using IdxRangeTor2 = IdxRange<GridTor2>;
using IdxRangeVpar = IdxRange<GridVpar>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeDim6 = IdxRange<GridDim6>;
using IdxRangeTor2D = IdxRange<GridTor1, GridTor2>;
using IdxRangeV2D = IdxRange<GridVpar, GridMu>;
using IdxRangeSpTor2DV2DDim6 = IdxRange<Species, GridTor1, GridTor2, GridVpar, GridMu, GridDim6>;

// Field type aliases
template <class ElementType>
using FieldMemSpTor2DV2DDim6 = FieldMem<ElementType, IdxRangeSpTor2DV2DDim6>;
using DFieldMemSpTor2DV2DDim6 = FieldMemSpTor2DV2DDim6<double>;

template <class ElementType>
using FieldSpTor2DV2DDim6 = Field<ElementType, IdxRangeSpTor2DV2DDim6>;
using DFieldSpTor2DV2DDim6 = FieldSpTor2DV2DDim6<double>;

template <class ElementType>
using FieldMemSpTor2DV2DDim6_host = host_t<FieldMem<ElementType, IdxRangeSpTor2DV2DDim6>>;
using DFieldMemSpTor2DV2DDim6_host = FieldMemSpTor2DV2DDim6_host<double>;

