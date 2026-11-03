// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/kernels/splines.hpp>

#include "geometry_vpar_mu.hpp"

// Splines definition
int constexpr BSDegreeVpar = 3;
int constexpr BSDegreeMu = 3;

bool constexpr BsplineOnUniformCellsVpar = true;
bool constexpr BsplineOnUniformCellsMu = true;

struct BSplinesVpar
    : std::conditional_t<
              BsplineOnUniformCellsVpar,
              ddc::UniformBSplines<Vpar, BSDegreeVpar, Vpar::PERIODIC>,
              ddc::NonUniformBSplines<Vpar, BSDegreeVpar, Vpar::PERIODIC>>
{
};
struct BSplinesMu
    : std::conditional_t<
              BsplineOnUniformCellsMu,
              ddc::UniformBSplines<Mu, BSDegreeMu, Mu::PERIODIC>,
              ddc::NonUniformBSplines<Mu, BSDegreeMu, Mu::PERIODIC>>
{
};
ddc::SplineBuilderClosure constexpr SplineVparClosure = ddc::SplineBuilderClosure::HERMITE;
ddc::SplineBuilderClosure constexpr SplineMuClosure = ddc::SplineBuilderClosure::HERMITE;

using SplineInterpPointsVpar
        = ddc::GrevilleInterpolationPoints<BSplinesVpar, SplineVparClosure, SplineVparClosure>;
using SplineInterpPointsMu
        = ddc::GrevilleInterpolationPoints<BSplinesMu, SplineMuClosure, SplineMuClosure>;

using SplineVparBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVpar,
        GridVpar,
        SplineVparClosure,
        SplineVparClosure,
        ddc::SplineSolver::LAPACK>;
using SplineVparEvaluator = ddc::SplineEvaluatorND<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        ddc::detail::TypeSeq<BSplinesVpar>,
        ddc::detail::TypeSeq<GridVpar>,
        ddc::detail::TypeSeq<
                ddc::ConstantExtrapolationRule<Vpar>,
                ddc::ConstantExtrapolationRule<Vpar>>>;

using SplineMuBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesMu,
        GridMu,
        SplineMuClosure,
        SplineMuClosure,
        ddc::SplineSolver::LAPACK>;
using SplineMuEvaluator = ddc::SplineEvaluatorND<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        ddc::detail::TypeSeq<BSplinesMu>,
        ddc::detail::TypeSeq<GridMu>,
        ddc::detail::
                TypeSeq<ddc::ConstantExtrapolationRule<Mu>, ddc::ConstantExtrapolationRule<Mu>>>;
