#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "ddc_aliases.hpp"

/*
 * @file geometry.hpp
 *
 */

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
using CoordVpar = Coord<Vpar>;
using CoordMu = Coord<Mu>;

// Splines definition
int constexpr BSDegreeVpar = 3;
int constexpr BSDegreeMu = 3;

bool constexpr BsplineOnUniformCellsVpar = true;
bool constexpr BsplineOnUniformCellsMu = true;

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
ddc::BoundCond constexpr SplineVparBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineMuBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsVpar
        = ddc::GrevilleInterpolationPoints<BSplinesVpar, SplineVparBoundary, SplineVparBoundary>;
using SplineInterpPointsMu
        = ddc::GrevilleInterpolationPoints<BSplinesMu, SplineMuBoundary, SplineMuBoundary>;

struct GridVpar : SplineInterpPointsVpar::interpolation_discrete_dimension_type
{
};
struct GridMu : SplineInterpPointsMu::interpolation_discrete_dimension_type
{
};

// Idx = index of the point in the point sampling
using IdxVpar = Idx<GridVpar>;
using IdxMu = Idx<GridMu>;
using IdxVparMu = Idx<GridVpar, GridMu>;
using IdxSpVparMu = Idx<Species, GridVpar, GridMu>;

// IdxStep = number of grid points between points in a sampling
using IdxStepVpar = IdxStep<GridVpar>;
using IdxStepMu = IdxStep<GridMu>;
using IdxStepVparMu = IdxStep<GridVpar, GridMu>;
using IdxStepSpVparMu = IdxStep<Species, GridVpar, GridMu>;

// IdxRange = to describe the wole index range (or a sub-index range)
using IdxRangeVpar = IdxRange<GridVpar>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeVparMu = IdxRange<GridVpar, GridMu>;
using IdxRangeSpVparMu = IdxRange<Species, GridVpar, GridMu>;

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


/* A OTER
New notation:
- RDimX --> X
- IDimX --> GridX
- IndexX --> IdxX
- IVectX --> IdxStepX
- IDomX --> IdxRangeX
- DBSsplineX --> SplineCoeffX
*/
