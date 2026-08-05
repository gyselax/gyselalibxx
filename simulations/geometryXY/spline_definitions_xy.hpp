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

ddc::SplineBuilderClosure constexpr SplineXClosure = ddc::SplineBuilderClosure::PERIODIC;
ddc::SplineBuilderClosure constexpr SplineYClosure = ddc::SplineBuilderClosure::PERIODIC;

using SplineXExtrapolation = ExtrapolationRule::Periodic;
using SplineYExtrapolation = ExtrapolationRule::Periodic;

// IDim initialisers
using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXClosure, SplineXClosure>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYClosure, SplineYClosure>;


// SplineBuilder and SplineEvaluator definitions
using SplineXInterpolator = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesX,
        GridX,
        SplineXExtrapolation,
        SplineBoundaryClosures<SplineXClosure, SplineXClosure>>;

using SplineYInterpolator = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesY,
        GridY,
        SplineYExtrapolation,
        SplineBoundaryClosures<SplineYClosure, SplineYClosure>>;

// Spline index range
using IdxRangeBSX = IdxRange<BSplinesX>;
using IdxRangeBSY = IdxRange<BSplinesY>;
using IdxRangeBSXY = IdxRange<BSplinesX, BSplinesY>;

template <class ElementType>
using BSConstFieldXY = Field<ElementType const, IdxRangeBSXY>;
using DBSConstFieldXY = BSConstFieldXY<double>;
