#pragma once
#include <ddc/kernels/splines.hpp>

#include "geometry.hpp"

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
auto constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;

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

template <class ElementType>
using BSFieldMemX = FieldMem<ElementType, IdxRangeBSX>;

using DBSFieldMemX = BSFieldMemX<double>;

template <class ElementType>
using BSFieldX = Field<ElementType, IdxRangeBSX>;

using DBSFieldX = BSFieldX<double>;

template <class ElementType>
using BSConstFieldX = ConstField<ElementType, IdxRangeBSX>;

using DBSConstFieldX = BSConstFieldX<double>;
