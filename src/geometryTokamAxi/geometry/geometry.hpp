// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

/*
 * @file geometry.hpp
 *
 */

/**
 * @brief Define non periodic real R dimension.
 */
struct R
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define periodic real Theta dimension.
 */
struct Theta
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, periodic.
     */
    static bool constexpr PERIODIC = true;
};
/**
 * @brief Define non periodic parallel velocity @f$v_\parallel@f$.
 */
struct Vpar
{
    /**
     * The periodicity of the parallel velocity.
     * This is a compile time constant which allows the code to avoid
     * unnecessary if conditions.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Define non periodic magnetic momentum @f$\mu@f$.
 */
struct Mu
{
    /**
     * The periodicity of the magnetic momentum.
     * This is a compile time constant which allows the code to avoid
     * unnecessary if conditions.
     */
    static bool constexpr PERIODIC = false;
};

// Coord = position of a coordinate in the vector space
using CoordR = Coord<R>;
using CoordTheta = Coord<Theta>;
using CoordVpar = Coord<Vpar>;
using CoordMu = Coord<Mu>;

// Splines definition
int constexpr BSDegreeR = 3;
int constexpr BSDegreeTheta = 3;
int constexpr BSDegreeVpar = 3;
int constexpr BSDegreeMu = 3;

bool constexpr BsplineOnUniformCellsR = false;
bool constexpr BsplineOnUniformCellsTheta = true;
bool constexpr BsplineOnUniformCellsVpar = true;
bool constexpr BsplineOnUniformCellsMu = true;

struct BSplinesR
    : std::conditional_t<
              BsplineOnUniformCellsR,
              ddc::UniformBSplines<R, BSDegreeR>,
              ddc::NonUniformBSplines<R, BSDegreeR>>
{
};
struct BSplinesTheta
    : std::conditional_t<
              BsplineOnUniformCellsTheta,
              ddc::UniformBSplines<Theta, BSDegreeTheta>,
              ddc::NonUniformBSplines<Theta, BSDegreeTheta>>
{
};
struct BSplinesVpar
    : std::conditional_t<
              BsplineOnUniformCellsVpar,
              ddc::UniformBSplines<Vpar, BSDegreeVpar>,
              ddc::NonUniformBSplines<Vpar, BSDegreeVpar>>
{
};
struct BSplinesMu
    : std::conditional_t<
              BsplineOnUniformCellsMu,
              ddc::UniformBSplines<Mu, BSDegreeMu>,
              ddc::NonUniformBSplines<Mu, BSDegreeMu>>
{
};
ddc::BoundCond constexpr SplineRBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineThetaBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVparBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineMuBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaBoundary, SplineThetaBoundary>;
using SplineInterpPointsVpar
        = ddc::GrevilleInterpolationPoints<BSplinesVpar, SplineVparBoundary, SplineVparBoundary>;
using SplineInterpPointsMu
        = ddc::GrevilleInterpolationPoints<BSplinesMu, SplineMuBoundary, SplineMuBoundary>;

struct GridR : SplineInterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : SplineInterpPointsTheta::interpolation_discrete_dimension_type
{
};
struct GridVpar : SplineInterpPointsVpar::interpolation_discrete_dimension_type
{
};
struct GridMu : SplineInterpPointsMu::interpolation_discrete_dimension_type
{
};

// Idx = index of the point in the point sampling
using IdxR = Idx<GridR>;
using IdxTheta = Idx<GridTheta>;
using IdxVpar = Idx<GridVpar>;
using IdxCS = Idx<GridR, GridTheta>;
using IdxMu = Idx<GridMu>;
using IdxVparMu = Idx<GridVpar, GridMu>;
using IdxSpVparMu = Idx<Species, GridVpar, GridMu>;

// IdxStep = number of grid points between points in a sampling
using IdxStepR = IdxStep<GridR>;
using IdxStepTheta = IdxStep<GridTheta>;
using IdxStepCS = IdxStep<GridCS>;
using IdxStepVpar = IdxStep<GridVpar>;
using IdxStepMu = IdxStep<GridMu>;
using IdxStepVparMu = IdxStep<GridVpar, GridMu>;
using IdxStepSpVparMu = IdxStep<Species, GridVpar, GridMu>;

// IdxRange = to describe the wole index range (or a sub-index range)
using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeCS = IdxRange<GridCS>;
using IdxRangeVpar = IdxRange<GridVpar>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeVparMu = IdxRange<GridVpar, GridMu>;
using IdxRangeSpVparMu = IdxRange<Species, GridVpar, GridMu>;

// template for the fields
// --> For FieldMem template
template <class ElementType>
using FieldMemR = FieldMem<ElementType, IdxRangeR>;
using DFieldMemR = FieldMemR<double>;

template <class ElementType>
using FieldMemTheta = FieldMem<ElementType, IdxRangeTheta>;
using DFieldMemTheta = FieldMemTheta<double>;

template <class ElementType>
using FieldMemCS = FieldMem<ElementType, IdxRangeCS>;
using DFieldMemCS = FieldMemCS<double>;

template <class ElementType>
using FieldMemVpar = FieldMem<ElementType, IdxRangeVpar>;
using DFieldMemVpar = FieldMemVpar<double>;

template <class ElementType>
using FieldMemMu = FieldMem<ElementType, IdxRangeMu>;
using DFieldMemMu = FieldMemMu<double>;

template <class ElementType>
using FieldMemVparMu = FieldMem<ElementType, IdxRangeVparMu>;
using DFieldMemVparMu = FieldMemVparMu<double>;

// --> For Field template
template <class ElementType>
using FieldMemSpVparMu = FieldMem<ElementType, IdxRangeSpVparMu>;
using DFieldMemSpVparMu = FieldMemSpVparMu<double>;

template <class ElementType>
using FieldR = Field<ElementType, IdxRangeR>;
using DFieldR = FieldR<double>;

template <class ElementType>
using FieldTheta = Field<ElementType, IdxRangeTheta>;
using DFieldTheta = FieldTheta<double>;

template <class ElementType>
using FieldCS = Field<ElementType, IdxRangeCS>;
using DFieldCS = FieldCS<double>;

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

// --> For ConstField template
template <class ElementType>
using ConstFieldR = ConstField<ElementType, IdxRangeR>;
using DConstFieldR = ConstFieldR<double>;

template <class ElementType>
using ConstFieldTheta = ConstField<ElementType, IdxRangeTheta>;
using DConstFieldTheta = ConstFieldTheta<double>;

template <class ElementType>
using ConstFieldCS = ConstField<ElementType, IdxRangeCS>;
using DConstFieldCS = ConstFieldCS<double>;

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
