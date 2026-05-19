// SPDX-License-Identifier: MIT

#pragma once

#include "geometry_xyvxvy.hpp"
#include "lagrange_interpolation.hpp"


int constexpr LDegreeX = 3;
int constexpr LDegreeY = 3;

int constexpr LDegreeVx = 3;
int constexpr LDegreeVy = 3;

bool constexpr LagrangeOnUniformCellsX = true;
bool constexpr LagrangeOnUniformCellsY = true;

bool constexpr LagrangeOnUniformCellsVx = true;
bool constexpr LagrangeOnUniformCellsVy = true;

struct LagrangeX
    : std::conditional_t<
              LagrangeOnUniformCellsX,
              UniformLagrangeBasis<X, LDegreeX>,
              NonUniformLagrangeBasis<X, LDegreeX>>
{
};

struct LagrangeY
    : std::conditional_t<
              LagrangeOnUniformCellsY,
              UniformLagrangeBasis<Y, LDegreeY>,
              NonUniformLagrangeBasis<Y, LDegreeY>>
{
};

struct LagrangeVx
    : std::conditional_t<
              LagrangeOnUniformCellsVx,
              UniformLagrangeBasis<Vx, LDegreeVx>,
              NonUniformLagrangeBasis<Vx, LDegreeVx>>
{
};

struct LagrangeVy
    : std::conditional_t<
              LagrangeOnUniformCellsVy,
              UniformLagrangeBasis<Vy, LDegreeVy>,
              NonUniformLagrangeBasis<Vy, LDegreeVy>>
{
};

ddc::BoundCond constexpr LagrangeYBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr LagrangeVyBoundary = ddc::BoundCond::HOMOGENEOUS_HERMITE;

// SplineBuilder and SplineEvaluator definition
using LagrangeInterpolatorX
        = LagrangeInterpolator<Kokkos::DefaultExecutionSpace, LagrangeX, GridX, PERIODIC, PERIODIC>;
using LagrangeInterpolatorY
        = LagrangeInterpolator<Kokkos::DefaultExecutionSpace, LagrangeY, GridY, PERIODIC, PERIODIC>;

using LagrangeInterpolatorVx = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagrangeVx,
        GridVx,
        CONSTANT,
        CONSTANT>;
using LagrangeInterpolatorVy = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagrangeVy,
        GridVy,
        CONSTANT,
        CONSTANT>;

using IdxRangeLY = IdxRange<LagrangeY>;
using IdxRangeLXY = IdxRange<LagrangeX, LagrangeY>;
using IdxRangeLVy = IdxRange<LagrangeVy>;
using IdxRangeLVxVy = IdxRange<LagrangeVx, LagrangeVy>;

template <class ElementType>
using LConstFieldXY = Field<ElementType const, IdxRangeLXY>;
using DLConstFieldXY = LConstFieldXY<double>;
