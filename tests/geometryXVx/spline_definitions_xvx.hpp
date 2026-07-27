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

auto constexpr SplineXClosure
        = X::PERIODIC ? ddc::SplineBuilderClosure::PERIODIC : ddc::SplineBuilderClosure::GREVILLE;
auto constexpr SplineVxClosure = ddc::SplineBuilderClosure::HERMITE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXClosure, SplineXClosure>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxClosure, SplineVxClosure>;

using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXClosure,
        SplineXClosure,
        ddc::SplineSolver::LAPACK>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        std::conditional_t<
                X::PERIODIC,
                ddc::PeriodicExtrapolationRule<X>,
                ddc::ConstantExtrapolationRule<X>>,
        std::conditional_t<
                X::PERIODIC,
                ddc::PeriodicExtrapolationRule<X>,
                ddc::ConstantExtrapolationRule<X>>>;
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        SplineVxClosure,
        SplineVxClosure,
        ddc::SplineSolver::LAPACK>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        ddc::ConstantExtrapolationRule<Vx>,
        ddc::ConstantExtrapolationRule<Vx>>;

using XExtrapRule
        = std::conditional_t<X::PERIODIC, ExtrapolationRule::Periodic, ExtrapolationRule::Constant>;

using SplineInterpolatorX = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesX,
        GridX,
        ddc::detail::TypeSeq<XExtrapRule, XExtrapRule>,
        SplineXClosure,
        SplineXClosure>;

using SplineInterpolatorVx = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesVx,
        GridVx,
        ddc::detail::TypeSeq<ExtrapolationRule::Constant, ExtrapolationRule::Constant>,
        SplineVxClosure,
        SplineVxClosure>;

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
