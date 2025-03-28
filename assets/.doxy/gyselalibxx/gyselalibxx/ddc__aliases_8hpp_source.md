

# File ddc\_aliases.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**ddc\_aliases.hpp**](ddc__aliases_8hpp.md)

[Go to the documentation of this file](ddc__aliases_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

template <class... Dims>
using Coord = ddc::Coordinate<Dims...>;

template <class... GridTypes>
using Idx = ddc::DiscreteElement<GridTypes...>;

template <class... GridTypes>
using IdxStep = ddc::DiscreteVector<GridTypes...>;

template <class... GridTypes>
using IdxRange = ddc::DiscreteDomain<GridTypes...>;

template <
        class ElementType,
        class IdxRange,
        class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
using FieldMem = ddc::Chunk<ElementType, IdxRange, ddc::KokkosAllocator<ElementType, MemSpace>>;

template <class IdxRange, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
using DFieldMem = FieldMem<double, IdxRange, MemSpace>;

template <
        class ElementType,
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using Field = ddc::ChunkSpan<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;

template <
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DField = Field<double, IdxRange, MemorySpace, LayoutStridedPolicy>;

template <
        class ElementType,
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using ConstField = ddc::ChunkView<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;

template <
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DConstField = ConstField<double, IdxRange, MemorySpace, LayoutStridedPolicy>;

template <class Dim>
using UniformGridBase = ddc::UniformPointSampling<Dim>;

template <class Dim>
using NonUniformGridBase = ddc::NonUniformPointSampling<Dim>;


//----------------------------------------------
//  Class for batched 1D spline builder
//----------------------------------------------
namespace detail {
template <
        class ExecSpace,
        class BSplines,
        class InterpolationGrid,
        ddc::BoundCond BCLBound,
        ddc::BoundCond BCUBound,
        ddc::SplineSolver SolverType,
        class IdxRangeType>
class GetSplineBatchedBuilder1D;


template <
        class ExecSpace,
        class BSplines,
        class InterpolationGrid,
        ddc::BoundCond BCLBound,
        ddc::BoundCond BCUBound,
        ddc::SplineSolver SolverType,
        class... Grid1D>
class GetSplineBatchedBuilder1D<
        ExecSpace,
        BSplines,
        InterpolationGrid,
        BCLBound,
        BCUBound,
        SolverType,
        IdxRange<Grid1D...>>
{
public:
    using type = ddc::SplineBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplines,
            InterpolationGrid,
            BCLBound,
            BCUBound,
            SolverType,
            Grid1D...>;
};


//----------------------------------------------
//  Class for a batched 2D splines builder
//----------------------------------------------
template <
        class ExecSpace,
        class BSplinesX1,
        class BSplinesX2,
        class InterpolationGridX1,
        class InterpolationGridX2,
        ddc::BoundCond BCLBoundX1,
        ddc::BoundCond BCUBoundX1,
        ddc::BoundCond BCLBoundX2,
        ddc::BoundCond BCUBoundX2,
        ddc::SplineSolver SolverType,
        class IdxRangeType>
class GetSplineBatchedBuilder2D;

template <
        class ExecSpace,
        class BSplinesX1,
        class BSplinesX2,
        class InterpolationGridX1,
        class InterpolationGridX2,
        ddc::BoundCond BCLBoundX1,
        ddc::BoundCond BCUBoundX1,
        ddc::BoundCond BCLBoundX2,
        ddc::BoundCond BCUBoundX2,
        ddc::SplineSolver SolverType,
        class... Grid1D>
class GetSplineBatchedBuilder2D<
        ExecSpace,
        BSplinesX1,
        BSplinesX2,
        InterpolationGridX1,
        InterpolationGridX2,
        BCLBoundX1,
        BCUBoundX1,
        BCLBoundX2,
        BCUBoundX2,
        SolverType,
        IdxRange<Grid1D...>>
{
public:
    using type = ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplinesX1,
            BSplinesX2,
            InterpolationGridX1,
            InterpolationGridX2,
            BCLBoundX1,
            BCUBoundX1,
            BCLBoundX2,
            BCUBoundX2,
            SolverType,
            Grid1D...>;
};


//----------------------------------------------
//  Class for batched 1D splines evaluator
//----------------------------------------------
template <
        class ExecSpace,
        class BSplines,
        class InterpolationGrid,
        class LowerExtrapolationRule,
        class UpperExtrapolationRule,
        class IdxRangeType>
class GetSplineBatchedEvaluator1D;

template <
        class ExecSpace,
        class BSplines,
        class InterpolationGrid,
        class LowerExtrapolationRule,
        class UpperExtrapolationRule,
        class... Grid1D>
class GetSplineBatchedEvaluator1D<
        ExecSpace,
        BSplines,
        InterpolationGrid,
        LowerExtrapolationRule,
        UpperExtrapolationRule,
        IdxRange<Grid1D...>>
{
public:
    using type = ddc::SplineEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplines,
            InterpolationGrid,
            LowerExtrapolationRule,
            UpperExtrapolationRule,
            Grid1D...>;
};


//----------------------------------------------
//  Class for 2D splines evaluator
//----------------------------------------------
template <
        class ExecSpace,
        class BSplinesX1,
        class BSplinesX2,
        class InterpolationGridX1,
        class InterpolationGridX2,
        class LowerExtrapolationRuleX1,
        class UpperExtrapolationRuleX1,
        class LowerExtrapolationRuleX2,
        class UpperExtrapolationRuleX2,
        class IdxRangeType>
class GetSplineBatchedEvaluator2D;

template <
        class ExecSpace,
        class BSplinesX1,
        class BSplinesX2,
        class InterpolationGridX1,
        class InterpolationGridX2,
        class LowerExtrapolationRuleX1,
        class UpperExtrapolationRuleX1,
        class LowerExtrapolationRuleX2,
        class UpperExtrapolationRuleX2,
        class... Grid1D>
class GetSplineBatchedEvaluator2D<
        ExecSpace,
        BSplinesX1,
        BSplinesX2,
        InterpolationGridX1,
        InterpolationGridX2,
        LowerExtrapolationRuleX1,
        UpperExtrapolationRuleX1,
        LowerExtrapolationRuleX2,
        UpperExtrapolationRuleX2,
        IdxRange<Grid1D...>>
{
public:
    using type = ddc::SplineEvaluator2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            BSplinesX1,
            BSplinesX2,
            InterpolationGridX1,
            InterpolationGridX2,
            LowerExtrapolationRuleX1,
            UpperExtrapolationRuleX1,
            LowerExtrapolationRuleX2,
            UpperExtrapolationRuleX2,
            Grid1D...>;
};


} // namespace detail


//----------------------------------------------------
//  Aliases for batched  spline builder and evaluator
//----------------------------------------------------
template <
        class ExecSpace,
        class BSplines,
        class InterpolationGrid,
        ddc::BoundCond BCLBound,
        ddc::BoundCond BCUBound,
        ddc::SplineSolver SolverType,
        class IdxRangeType>
using get_spline_batched_builder1d_t = typename detail::GetSplineBatchedBuilder1D<
        ExecSpace,
        BSplines,
        InterpolationGrid,
        BCLBound,
        BCUBound,
        SolverType,
        IdxRangeType>::type;

template <
        class ExecSpace,
        class BSplinesX1,
        class BSplinesX2,
        class InterpolationGridX1,
        class InterpolationGridX2,
        ddc::BoundCond BCLBoundX1,
        ddc::BoundCond BCUBoundX1,
        ddc::BoundCond BCLBoundX2,
        ddc::BoundCond BCUBoundX2,
        ddc::SplineSolver SolverType,
        class IdxRangeType>
using get_spline_batched_builder2d_t = typename detail::GetSplineBatchedBuilder2D<
        ExecSpace,
        BSplinesX1,
        BSplinesX2,
        InterpolationGridX1,
        InterpolationGridX2,
        BCLBoundX1,
        BCUBoundX1,
        BCLBoundX2,
        BCUBoundX2,
        SolverType,
        IdxRangeType>::type;

template <
        class ExecSpace,
        class BSplines,
        class InterpolationGrid,
        class LowerExtrapolationRule,
        class UpperExtrapolationRule,
        class IdxRangeType>
using get_spline_batched_evaluator1d_t = typename detail::GetSplineBatchedEvaluator1D<
        ExecSpace,
        BSplines,
        InterpolationGrid,
        LowerExtrapolationRule,
        UpperExtrapolationRule,
        IdxRangeType>::type;

template <
        class ExecSpace,
        class BSplinesX1,
        class BSplinesX2,
        class InterpolationGridX1,
        class InterpolationGridX2,
        class LowerExtrapolationRuleX1,
        class UpperExtrapolationRuleX1,
        class LowerExtrapolationRuleX2,
        class UpperExtrapolationRuleX2,
        class IdxRangeType>
using get_spline_batched_evaluator2d_t = typename detail::GetSplineBatchedEvaluator2D<
        ExecSpace,
        BSplinesX1,
        BSplinesX2,
        InterpolationGridX1,
        InterpolationGridX2,
        LowerExtrapolationRuleX1,
        UpperExtrapolationRuleX1,
        LowerExtrapolationRuleX2,
        UpperExtrapolationRuleX2,
        IdxRangeType>::type;
```


