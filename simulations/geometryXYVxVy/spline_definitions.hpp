// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/kernels/splines.hpp>

#include "geometry.hpp"


int constexpr BSDegreeX = 3;
int constexpr BSDegreeY = 3;

int constexpr BSDegreeVx = 3;
int constexpr BSDegreeVy = 3;

bool constexpr BsplineOnUniformCellsX = true;
bool constexpr BsplineOnUniformCellsY = true;

bool constexpr BsplineOnUniformCellsVx = true;
bool constexpr BsplineOnUniformCellsVy = true;

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

struct BSplinesVx
    : std::conditional_t<
              BsplineOnUniformCellsVx,
              ddc::UniformBSplines<Vx, BSDegreeVx>,
              ddc::NonUniformBSplines<Vx, BSDegreeVx>>
{
};
struct BSplinesVy
    : std::conditional_t<
              BsplineOnUniformCellsVy,
              ddc::UniformBSplines<Vy, BSDegreeVy>,
              ddc::NonUniformBSplines<Vy, BSDegreeVy>>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;
ddc::BoundCond constexpr SplineVyBoundary = ddc::BoundCond::HERMITE;

// IDim initialisers
using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;
using SplineInterpPointsVy
        = ddc::GrevilleInterpolationPoints<BSplinesVy, SplineVyBoundary, SplineVyBoundary>;

// SplineBuilder and SplineEvaluator definition
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>>;
using SplineYBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineYEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        ddc::PeriodicExtrapolationRule<Y>,
        ddc::PeriodicExtrapolationRule<Y>>;
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        ddc::ConstantExtrapolationRule<Vx>,
        ddc::ConstantExtrapolationRule<Vx>>;
using SplineVyBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVy,
        GridVy,
        SplineVyBoundary,
        SplineVyBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineVyEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVy,
        GridVy,
        ddc::ConstantExtrapolationRule<Vy>,
        ddc::ConstantExtrapolationRule<Vy>>;

using IdxRangeBSX = IdxRange<BSplinesX>;
using IdxRangeBSY = IdxRange<BSplinesY>;
using IdxRangeBSXY = IdxRange<BSplinesX, BSplinesY>;
using IdxRangeBSVx = IdxRange<BSplinesVx>;
using IdxRangeBSVy = IdxRange<BSplinesVy>;
using IdxRangeBSVxVy = IdxRange<BSplinesVx, BSplinesVy>;

template <class ElementType>
using BSConstFieldXY = Field<ElementType const, IdxRangeBSXY>;
using DBSConstFieldXY = BSConstFieldXY<double>;
