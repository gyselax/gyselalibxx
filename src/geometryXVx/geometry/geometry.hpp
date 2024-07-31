// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>
#include <ddc/kernels/splines.hpp>

#include <ddc_helper.hpp>
#include <moments.hpp>
#include <species_info.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief A class which describes the real space in the spatial X direction
 */
struct X
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
#ifdef PERIODIC_RDIMX
    static bool constexpr PERIODIC = true;
#else
    static bool constexpr PERIODIC = false;
#endif
};

/**
 * @brief A class which describes the real space in the X-velocity direction
 */
struct Vx
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};

/**
 * @brief A class which describes the real space in the temporal direction
 */
struct T
{
    /**
     * @brief A boolean indicating if the dimension is periodic.
     */
    static bool constexpr PERIODIC = false;
};


using CoordT = Coord<T>;
using CoordX = Coord<X>;

using CoordVx = Coord<Vx>;

using CoordXVx = Coord<X, Vx>;

int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsVx = true;

struct BSplinesX
    : std::conditional_t<
              BsplineOnUniformCellsX,
              ddc::UniformBSplines<X, BSDegreeX>,
              ddc::NonUniformBSplines<X, BSDegreeX>>
{
};
struct BSplinesVx
    : std::conditional_t<
              BsplineOnUniformCellsVx,
              ddc::UniformBSplines<Vx, BSDegreeVx>,
              ddc::NonUniformBSplines<Vx, BSDegreeVx>>
{
};

auto constexpr SplineXBoundary = X::PERIODIC ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE;
auto constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridVx : SplineInterpPointsVx::interpolation_discrete_dimension_type
{
};

using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridVx>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
#ifdef PERIODIC_RDIMX
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>,
#else
        ddc::ConstantExtrapolationRule<X>,
        ddc::ConstantExtrapolationRule<X>,
#endif
        GridX,
        GridVx>;
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridVx>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        ddc::ConstantExtrapolationRule<Vx>,
        ddc::ConstantExtrapolationRule<Vx>,
        GridX,
        GridVx>;
using SplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX>;
using SplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
#ifdef PERIODIC_RDIMX
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>,
#else
        ddc::ConstantExtrapolationRule<X>,
        ddc::ConstantExtrapolationRule<X>,
#endif
        GridX>;
using SplineVxBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::LAPACK,
        GridVx>;
using SplineVxEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        ddc::ConstantExtrapolationRule<Vx>,
        ddc::ConstantExtrapolationRule<Vx>,
        GridVx>;

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



using BSIdxRangeX = IdxRange<BSplinesX>;

using BSIdxRangeVx = IdxRange<BSplinesVx>;



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
using BSFieldMemX = FieldMem<ElementType, BSIdxRangeX>;


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

using DBSFieldMemX = BSFieldMemX<double>;


using DFieldMemSpMom = FieldMemSpMom<double>;

using DFieldMemSpMomX = FieldMemSpMomX<double>;

using DFieldMemSpVx = FieldMemSpVx<double>;

using DFieldMemSpX = FieldMemSpX<double>;

using DFieldMemSpXVx = FieldMemSpXVx<double>;



template <class ElementType>
using BSFieldX = Field<ElementType, BSIdxRangeX>;

template <class ElementType>
using FieldX = Field<ElementType, IdxRangeX>;

template <class ElementType>
using FieldVx = Field<ElementType, IdxRangeVx>;

template <class ElementType>
using FieldSpMomX = Field<ElementType, IdxRangeSpMomX>;

template <class ElementType>
using FieldSpVx = Field<ElementType, IdxRangeSpVx>;

template <class ElementType>
using FieldSpX = Field<ElementType, IdxRangeSpX>;

template <class ElementType>
using FieldSpXVx = Field<ElementType, IdxRangeSpXVx>;


using DBSFieldX = BSFieldX<double>;

using DFieldVx = FieldVx<double>;

using DFieldX = FieldX<double>;

using DFieldSpMomX = FieldSpMomX<double>;

using DFieldSpVx = FieldSpVx<double>;

using DFieldSpX = FieldSpX<double>;

using DFieldSpXVx = FieldSpXVx<double>;


template <class ElementType>
using ConstFieldVx = Field<ElementType const, IdxRangeVx>;

template <class ElementType>
using ConstFieldX = Field<ElementType const, IdxRangeX>;


template <class ElementType>
using BSConstFieldX = Field<ElementType const, BSIdxRangeX>;

template <class ElementType>
using ConstFieldSpMom = Field<ElementType const, IdxRangeSpMom>;

template <class ElementType>
using ConstFieldSpMomX = Field<ElementType const, IdxRangeSpMomX>;

template <class ElementType>
using ConstFieldSpVx = Field<ElementType const, IdxRangeSpVx>;

template <class ElementType>
using ConstFieldSpX = Field<ElementType const, IdxRangeSpX>;

template <class ElementType>
using ConstFieldSpXVx = Field<ElementType const, IdxRangeSpXVx>;



using DConstFieldVx = ConstFieldVx<double>;

using DConstFieldX = ConstFieldX<double>;


using DBSConstFieldX = BSConstFieldX<double>;

using DConstFieldSpMom = ConstFieldSpMom<double>;

using DConstFieldSpMomX = ConstFieldSpMomX<double>;

using DConstFieldSpX = ConstFieldSpX<double>;

using DConstFieldSpVx = ConstFieldSpVx<double>;

using DConstFieldSpXVx = ConstFieldSpXVx<double>;


/**
 * @brief A class providing aliases for useful subindex ranges of the geometry. It is used as template parameter for generic dimensionality-agnostic operators such as advections.
 */
class GeometryXVx
{
public:
    /**
     * @brief A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type.
     */
    template <class T>
    using velocity_dim_for = std::conditional_t<std::is_same_v<T, GridX>, GridVx, void>;

    /**
     * @brief A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type.
     */
    template <class T>
    using spatial_dim_for = std::conditional_t<std::is_same_v<T, GridVx>, GridX, void>;

    /**
     * @brief An alias for the spatial discrete index range type.
     */
    using SpatialIdxRange = IdxRangeX;

    /**
     * @brief An alias for the velocity discrete index range type.
     */
    using VelocityIdxRange = IdxRangeVx;


    // using FdistribuIdxRange = DiscreteDomain<DimSp, typename decltype(SpatialDDom), typename decltype(VelocityDDom)>(IdxRange());
    /**
     * @brief An alias for the whole distribution function discrete index range type.
     */
    using FdistribuIdxRange = IdxRangeSpXVx;
};
