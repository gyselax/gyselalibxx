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
              UniformLagrangeBasis<X, LDegreeX, Real>,
              NonUniformLagrangeBasis<X, LDegreeX, Real>>
{
};

struct LagrangeY
    : std::conditional_t<
              LagrangeOnUniformCellsY,
              UniformLagrangeBasis<Y, LDegreeY, Real>,
              NonUniformLagrangeBasis<Y, LDegreeY, Real>>
{
};

struct LagrangeVx
    : std::conditional_t<
              LagrangeOnUniformCellsVx,
              UniformLagrangeBasis<Vx, LDegreeVx, Real>,
              NonUniformLagrangeBasis<Vx, LDegreeVx, Real>>
{
};

struct LagrangeVy
    : std::conditional_t<
              LagrangeOnUniformCellsVy,
              UniformLagrangeBasis<Vy, LDegreeVy, Real>,
              NonUniformLagrangeBasis<Vy, LDegreeVy, Real>>
{
};

ddc::SplineBuilderClosure constexpr LagrangeYBoundary = ddc::SplineBuilderClosure::PERIODIC;
ddc::SplineBuilderClosure constexpr LagrangeVyBoundary
        = ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE;

// SplineBuilder and SplineEvaluator definition
using LagrangeInterpolatorX = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagrangeX,
        GridX,
        PERIODIC,
        PERIODIC,
        Real>;
using LagrangeInterpolatorY = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagrangeY,
        GridY,
        PERIODIC,
        PERIODIC,
        Real>;

using LagrangeInterpolatorVx = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagrangeVx,
        GridVx,
        CONSTANT,
        CONSTANT,
        Real>;
using LagrangeInterpolatorVy = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagrangeVy,
        GridVy,
        CONSTANT,
        CONSTANT,
        Real>;

using IdxRangeLY = IdxRange<LagrangeY>;
using IdxRangeLXY = IdxRange<LagrangeX, LagrangeY>;
using IdxRangeLVy = IdxRange<LagrangeVy>;
using IdxRangeLVxVy = IdxRange<LagrangeVx, LagrangeVy>;

template <class ElementType>
using LConstFieldXY = Field<ElementType const, IdxRangeLXY>;
