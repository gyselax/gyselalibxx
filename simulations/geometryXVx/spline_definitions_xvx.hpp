// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "geometry_xvx.hpp"
#include "spline_interpolation.hpp"

int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

#ifdef INPUT_MESH
bool constexpr BsplineOnUniformCellsX = false;
bool constexpr BsplineOnUniformCellsVx = false;
#else
bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsVx = true;
#endif

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
auto constexpr SplineVxBoundary = ddc::BoundCond::HOMOGENEOUS_HERMITE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;

ExtrapolationRule constexpr XExtrapRule = X::PERIODIC ? PERIODIC : CONSTANT;

using SplineInterpolatorX = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesX,
        GridX,
        XExtrapRule,
        XExtrapRule,
        SplineXBoundary,
        SplineXBoundary>;

using SplineInterpolatorVx = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesVx,
        GridVx,
        CONSTANT,
        CONSTANT,
        SplineVxBoundary,
        SplineVxBoundary>;

using IdxRangeBSX = IdxRange<BSplinesX>;

using IdxRangeBSVx = IdxRange<BSplinesVx>;

template <class ElementType>
using BSFieldMemX = FieldMem<ElementType, IdxRangeBSX>;

using DBSFieldMemX = BSFieldMemX<double>;

template <class ElementType>
using BSFieldX = Field<ElementType, IdxRangeBSX>;

using DBSFieldX = BSFieldX<double>;

template <class ElementType>
using BSConstFieldX = ConstField<ElementType, IdxRangeBSX>;

using DBSConstFieldX = BSConstFieldX<double>;
