#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <ddc_helper.hpp>
#include <species_info.hpp>

/*
 * @file geometry.hpp
 *
 */

using IdxSp = IndexSp;
using IdxRangeSp = IDomainSp;

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

// ddc::Coordinate = position of a coordinate in the vector space
using CoordVpar = ddc::Coordinate<Vpar>;
using CoordMu = ddc::Coordinate<Mu>;

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

// ddc::DiscreteElement = index of the point in the point sampling
using IdxVpar = ddc::DiscreteElement<GridVpar>;
using IdxMu = ddc::DiscreteElement<GridMu>;
using IdxVparMu = ddc::DiscreteElement<GridVpar, GridMu>;
using IdxSpVparMu = ddc::DiscreteElement<IDimSp, GridVpar, GridMu>;

// ddc::DiscreteVector = number of grid points between points in a sampling
using IdxStepVpar = ddc::DiscreteVector<GridVpar>;
using IdxStepMu = ddc::DiscreteVector<GridMu>;
using IdxStepVparMu = ddc::DiscreteVector<GridVpar, GridMu>;
using IdxStepSpVparMu = ddc::DiscreteVector<IDimSp, GridVpar, GridMu>;

// ddc::DiscreteDomain = to describe the wole domain (or a sub-domain)
using IdxRangeVpar = ddc::DiscreteDomain<GridVpar>;
using IdxRangeMu = ddc::DiscreteDomain<GridMu>;
using IdxRangeVparMu = ddc::DiscreteDomain<GridVpar, GridMu>;
using IdxRangeSpVparMu = ddc::DiscreteDomain<IDimSp, GridVpar, GridMu>;

// template for the fields
template <class ElementType>
using FieldVpar = device_t<ddc::Chunk<ElementType, IdxRangeVpar>>;
using DFieldVpar = FieldVpar<double>;

template <class ElementType>
using FieldMu = device_t<ddc::Chunk<ElementType, IdxRangeMu>>;
using DFieldMu = FieldMu<double>;

template <class ElementType>
using FieldVparMu = device_t<ddc::Chunk<ElementType, IdxRangeVparMu>>;
using DFieldVparMu = FieldVparMu<double>;

template <class ElementType>
using FieldSpVparMu = device_t<ddc::Chunk<ElementType, IdxRangeSpVparMu>>;
using DFieldSpVparMu = FieldSpVparMu<double>;

template <class ElementType>
using SpanVparMu = device_t<ddc::ChunkSpan<ElementType, IdxRangeVparMu>>;
using DSpanVparMu = SpanVparMu<double>;

template <class ElementType>
using SpanSpVparMu = device_t<ddc::ChunkSpan<ElementType, IdxRangeSpVparMu>>;
using DSpanSpVparMu = SpanSpVparMu<double>;

template <class ElementType>
using ViewVparMu = device_t<ddc::ChunkView<ElementType, IdxRangeVparMu>>;
using DViewVparMu = ViewVparMu<double>;

template <class ElementType>
using ViewSpVparMu = device_t<ddc::ChunkView<ElementType, IdxRangeSpVparMu>>;
using DViewSpVparMu = ViewSpVparMu<double>;


/* A OTER
New notation:
- RDimX --> X
- IDimX --> GridX
- IndexX --> IdxX
- IVectX --> IdxStepX
- IDomX --> IdxRangeX
- DBSsplineX --> SplineCoeffX
*/
