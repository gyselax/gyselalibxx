// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "geometry_xy.hpp"
#include "spline_interpolation.hpp"


int constexpr BSDegreeX = 3;
int constexpr BSDegreeY = 3;

bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsY = true;

struct BSplinesX
    : std::conditional_t<
              BsplineOnUniformCellsX,
              ddc::UniformBSplines<X, BSDegreeX>,
              ddc::NonUniformBSplines<X, BSDegreeX>>
{
};
struct BSplinesY
    : std::conditional_t<
              BsplineOnUniformCellsY,
              ddc::UniformBSplines<Y, BSDegreeY>,
              ddc::NonUniformBSplines<Y, BSDegreeY>>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::PERIODIC;

ExtrapolationRule constexpr SplineXExtrapolation = ExtrapolationRule::PERIODIC;
ExtrapolationRule constexpr SplineYExtrapolation = ExtrapolationRule::PERIODIC;

// IDim initialisers
using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;


// SplineBuilder and SplineEvaluator definitions
using SplineXInterpolator = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesX,
        GridX,
        SplineXExtrapolation,
        SplineXExtrapolation,
        SplineXBoundary,
        SplineXBoundary>;

using SplineYInterpolator = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesY,
        GridY,
        SplineYExtrapolation,
        SplineYExtrapolation,
        SplineYBoundary,
        SplineYBoundary>;

// Spline index range
using IdxRangeBSX = IdxRange<BSplinesX>;
using IdxRangeBSY = IdxRange<BSplinesY>;
using IdxRangeBSXY = IdxRange<BSplinesX, BSplinesY>;

template <class ElementType>
using BSConstFieldXY = Field<ElementType const, IdxRangeBSXY>;
using DBSConstFieldXY = BSConstFieldXY<double>;
