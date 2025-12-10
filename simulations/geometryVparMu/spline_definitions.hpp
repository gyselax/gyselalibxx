// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/kernels/splines.hpp>
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

using SplineVparBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVpar,
        GridVpar,
        SplineVparBoundary,
        SplineVparBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineVparEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVpar,
        GridVpar,
        ddc::ConstantExtrapolationRule<Vpar>,
        ddc::ConstantExtrapolationRule<Vpar>>;

using SplineMuBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesMu,
        GridMu,
        SplineMuBoundary,
        SplineMuBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineMuEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesMu,
        GridMu,
        ddc::ConstantExtrapolationRule<Mu>,
        ddc::ConstantExtrapolationRule<Mu>>;
