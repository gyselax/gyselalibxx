// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/kernels/splines.hpp>

#include "../geometryXVx/spline_definitions_xvx.hpp"

#include "geometry_xyvxvy.hpp"


int constexpr BSDegreeY = 3;

int constexpr BSDegreeVy = 3;

bool constexpr BsplineOnUniformCellsY = true;

bool constexpr BsplineOnUniformCellsVy = true;

struct BSplinesY
    : std::conditional_t<
              BsplineOnUniformCellsY,
              ddc::UniformBSplines<Y, BSDegreeY>,
              ddc::NonUniformBSplines<Y, BSDegreeY>>
{
};

struct BSplinesVy
    : std::conditional_t<
              BsplineOnUniformCellsVy,
              ddc::UniformBSplines<Vy, BSDegreeVy>,
              ddc::NonUniformBSplines<Vy, BSDegreeVy>>
{
};

ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVyBoundary = ddc::BoundCond::HOMOGENEOUS_HERMITE;

// IDim initialisers
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
using SplineInterpPointsVy
        = ddc::GrevilleInterpolationPoints<BSplinesVy, SplineVyBoundary, SplineVyBoundary>;

// SplineBuilder and SplineEvaluator definition
using SplineInterpolatorY = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesY,
        GridY,
        PERIODIC,
        PERIODIC,
        SplineYBoundary,
        SplineYBoundary>;

using SplineInterpolatorVy = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesVy,
        GridVy,
        CONSTANT,
        CONSTANT,
        SplineVyBoundary,
        SplineVyBoundary>;

using IdxRangeBSY = IdxRange<BSplinesY>;
using IdxRangeBSXY = IdxRange<BSplinesX, BSplinesY>;
using IdxRangeBSVy = IdxRange<BSplinesVy>;
using IdxRangeBSVxVy = IdxRange<BSplinesVx, BSplinesVy>;

template <class ElementType>
using BSConstFieldXY = Field<ElementType const, IdxRangeBSXY>;
using DBSConstFieldXY = BSConstFieldXY<double>;
