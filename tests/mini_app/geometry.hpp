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
 * to initialize and read a 5D distribution function (species, tor1, tor2, tor3, vpar, mu)
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


struct Tor3
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
};
using Phi = Tor3;

// Velocity space dimensions
struct Vpar
{
    static bool constexpr PERIODIC = false;
};

struct Mu
{
    static bool constexpr PERIODIC = false;
};

using CoordTor1 = Coord<Tor1>;
using CoordTor2 = Coord<Tor2>;
using CoordTor3 = Coord<Tor3>;
using CoordVpar = Coord<Vpar>;
using CoordMu = Coord<Mu>;


// Splines definition
int constexpr BSDegreeTor1 = 3;
int constexpr BSDegreeTor2 = 3;
int constexpr BSDegreeTor3 = 3;
int constexpr BSDegreeVpar = 3;
int constexpr BSDegreeMu = 3;

struct BSplinesTor1 : ddc::UniformBSplines<Tor1, BSDegreeTor1>
{
};
struct BSplinesTor2 : ddc::UniformBSplines<Tor2, BSDegreeTor2>
{
};
struct BSplinesTor3 : ddc::UniformBSplines<Tor3, BSDegreeTor3>
{
};
struct BSplinesVpar : ddc::UniformBSplines<Vpar, BSDegreeVpar>
{
};
struct BSplinesMu : ddc::UniformBSplines<Mu, BSDegreeMu>
{
};

ddc::BoundCond constexpr SplineTor1Boundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineTor2Boundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineTor3Boundary = ddc::BoundCond::PERIODIC; 
ddc::BoundCond constexpr SplineVparBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineMuBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsTor1
        = ddc::GrevilleInterpolationPoints<BSplinesTor1, SplineTor1Boundary, SplineTor1Boundary>;
using SplineInterpPointsTor2
        = ddc::GrevilleInterpolationPoints<BSplinesTor2, SplineTor2Boundary, SplineTor2Boundary>;
using SplineInterpPointsTor3
        = ddc::GrevilleInterpolationPoints<BSplinesTor3, SplineTor3Boundary, SplineTor3Boundary>;
using SplineInterpPointsVpar
        = ddc::GrevilleInterpolationPoints<BSplinesVpar, SplineVparBoundary, SplineVparBoundary>;
using SplineInterpPointsMu
        = ddc::GrevilleInterpolationPoints<BSplinesMu, SplineMuBoundary, SplineMuBoundary>;


struct GridTor1 : SplineInterpPointsTor1::interpolation_discrete_dimension_type{};
struct GridTor2 : SplineInterpPointsTor2::interpolation_discrete_dimension_type{};
struct GridTor3 : SplineInterpPointsTor3::interpolation_discrete_dimension_type{};
struct GridVpar : SplineInterpPointsVpar::interpolation_discrete_dimension_type{};
struct GridMu : SplineInterpPointsMu::interpolation_discrete_dimension_type{};

// IdxRange definitions
using IdxRangeTor1 = IdxRange<GridTor1>;
using IdxRangeTor2 = IdxRange<GridTor2>;
using IdxRangeTor3 = IdxRange<GridTor3>;
using IdxRangeVpar = IdxRange<GridVpar>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeSpGrid = IdxRange<Species, GridTor1, GridTor2, GridTor3, GridVpar, GridMu>;
using IdxRangeSpVparMu = IdxRange<Species, GridVpar, GridMu>;
using IdxSpGrid = Idx<Species, GridTor1, GridTor2, GridTor3, GridVpar, GridMu>;
using IdxSpVparMu = Idx<Species, GridVpar, GridMu>;


template <class ElementType>
using FieldMemSpGrid = FieldMem<ElementType, IdxRangeSpGrid>;
using DFieldMemSpGrid = FieldMemSpGrid<double>;

template <class ElementType>
using FieldMemSpVparMu = FieldMem<ElementType, IdxRangeSpVparMu>;
using DFieldMemSpVparMu = FieldMemSpVparMu<double>;
