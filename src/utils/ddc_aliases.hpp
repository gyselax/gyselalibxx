// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

/**
 * This file contains aliases for DDC. The documentation for DDC can be found at <https://ddc.mdls.fr/>.
 * The names used in the DDC project are not always intuitive for mathematicians/physicists therefore
 * Gysela has chosen to introduce this file to provide names that are hopefully more intuitive.
 * The documentation for these concepts can be found in the Gysela documentation at:
 * <https://gyselax.github.io/gyselalibxx/docs_DDC_in_gyselalibxx.html>
 */

/// An alias describing the type of a coordinate (e.g. a coordinate in phase-space (x, vx)).
template <class... Dims>
using Coord = ddc::Coordinate<Dims...>;

/// An alias describing the type of an index that is used to access the values of a field defined on a grid.
template <class... GridTypes>
using Idx = ddc::DiscreteElement<GridTypes...>;

/// An alias describing the type of a distance between two indexes.
template <class... GridTypes>
using IdxStep = ddc::DiscreteVector<GridTypes...>;

/// An alias describing the type of an index range describing the subsection of a grid on which a field is defined.
template <class... GridTypes>
using IdxRange = ddc::DiscreteDomain<GridTypes...>;

/// An alias describing the type of an object which will allocate memory for a field when it is created.
template <
        class ElementType,
        class IdxRange,
        class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
using FieldMem = ddc::Chunk<ElementType, IdxRange, ddc::KokkosAllocator<ElementType, MemSpace>>;

/// An alias describing the type of an object which will allocate memory for a field of doubles when it is created.
template <class IdxRange, class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
using DFieldMem = FieldMem<double, IdxRange, MemSpace>;

/// An alias describing the type of a field defined on a grid (e.g. the electric field defined on the grid @f${x_0, x_1, .., x_N}@f$)
template <
        class ElementType,
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using Field = ddc::ChunkSpan<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;

/// An alias describing the type of a field of doubles defined on a grid (e.g. the electric field defined on the grid @f${x_0, x_1, .., x_N}@f$)
template <
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DField = Field<double, IdxRange, MemorySpace, LayoutStridedPolicy>;

/// An alias describing the type of a constant field defined on a grid (e.g. the electric field defined on the grid @f${x_0, x_1, .., x_N}@f$)
template <
        class ElementType,
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using ConstField = ddc::ChunkView<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;

/// An alias describing the type of a constant field of doubles defined on a grid (e.g. the electric field defined on the grid @f${x_0, x_1, .., x_N}@f$)
template <
        class IdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DConstField = ConstField<double, IdxRange, MemorySpace, LayoutStridedPolicy>;

/// An alias describing the type from which a uniform grid must inherit.
template <class GridType>
using UniformGridBase = ddc::UniformPointSampling<GridType>;

/// An alias describing the type from which a non-uniform grid must inherit.
template <class GridType>
using NonUniformGridBase = ddc::NonUniformPointSampling<GridType>;

/// A type describing a vector e.g. (E_x, E_y)
template <class ElementType, class... Dims>
class Vector : public ddc::detail::TaggedVector<ElementType, Dims...>
{
public:
    /**
     * @brief A constructor for a Vector to build from individual elements:
     * E.g. Vector<double, X, Y>(1.0, 2.0)
     *
     * or to build from component vectors
     * E.g. Vector<double, X, Y>(vec_x, vec_y)
     *
     * @param[in] p The element of the Vector.
     */
    template <class... Params>
    KOKKOS_FUNCTION Vector(Params... p) : ddc::detail::TaggedVector<ElementType, Dims...>(p...)
    {
    }
};

/// A type describing a vector whose elements are doubles e.g. (E_x, E_y)
template <class... Dims>
using DVector = Vector<double, Dims...>;


//----------------------------------------------
//  Aliases for batched 1D spline builder
//----------------------------------------------
/// A type describing a batched 1D spline builder
template<class ExecSpace, class BSplines, class InterpolationGrid, ddc::BoundCond BCLBound, ddc::BoundCond BCUBound, ddc::SplineSolver SolverType, class IdxRangeType>
class GetSplineBatchedBuilder1D;

template<class ExecSpace, class BSplines, class InterpolationGrid, ddc::BoundCond BCLBound, ddc::BoundCond BCUBound, ddc::SplineSolver SolverType, class... Grid1D>
class GetSplineBatchedBuilder1D<ExecSpace, BSplines, InterpolationGrid, BCLBound, BCUBound, SolverType, IdxRange<Grid1D...>> {
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

/// An alias describing a spline builder
template<class ExecSpace, class BSplines, class InterpolationGrid, ddc::BoundCond BCLBound, ddc::BoundCond BCUBound, ddc::SplineSolver SolverType, class IdxRangeType>
using get_spline_batched_builder1d_t = typename GetSplineBatchedBuilder1D<ExecSpace, BSplines, InterpolationGrid, BCLBound, BCUBound, SolverType, IdxRangeType>::type;


//----------------------------------------------
//  Aliases for a batched 2D splines builder
//----------------------------------------------
/// A type describing a batched 2D spline builder
template<class ExecSpace, class BSplinesX1, class BSplinesX2,
        class InterpolationGridX1, class InterpolationGridX2,
        ddc::BoundCond BCLBoundX1, ddc::BoundCond BCUBoundX1,
        ddc::BoundCond BCLBoundX2, ddc::BoundCond BCUBoundX2,
        ddc::SplineSolver SolverType, class IdxRangeType>
class GetSplineBatchedBuilder2D;

template<class ExecSpace, class BSplinesX1, class BSplinesX2, 
        class InterpolationGridX1, class InterpolationGridX2,
        ddc::BoundCond BCLBoundX1, ddc::BoundCond BCUBoundX1, 
        ddc::BoundCond BCLBoundX2, ddc::BoundCond BCUBoundX2, 
        ddc::SplineSolver SolverType, class... Grid1D>
class GetSplineBatchedBuilder2D<ExecSpace, BSplinesX1, BSplinesX2,
        InterpolationGridX1, InterpolationGridX2,
        BCLBoundX1, BCUBoundX1,
        BCLBoundX2, BCUBoundX2,
        SolverType, IdxRange<Grid1D...>> {
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

/// An alias describing a batched 2D spline builder
template<class ExecSpace, class BSplinesX1, class BSplinesX2, 
        class InterpolationGridX1, class InterpolationGridX2,
        ddc::BoundCond BCLBoundX1, ddc::BoundCond BCUBoundX1, 
        ddc::BoundCond BCLBoundX2, ddc::BoundCond BCUBoundX2, 
        ddc::SplineSolver SolverType, class IdxRangeType>
using get_spline_batched_builder2d_t = typename GetSplineBatchedBuilder2D<ExecSpace, 
        BSplinesX1, BSplinesX2,
        InterpolationGridX1, InterpolationGridX2, 
        BCLBoundX1, BCUBoundX1, 
        BCLBoundX2, BCUBoundX2, 
        SolverType, IdxRangeType>::type;


//----------------------------------------------
//  Aliases for batched 1D splines evaluator
//----------------------------------------------
template<class ExecSpace, class BSplines, class InterpolationGrid, 
        class LowerExtrapolationRule,  
        class UpperExtrapolationRule,  
        class IdxRangeType>
class GetSplineBatchedEvaluator1D;

template<class ExecSpace, class BSplines, class InterpolationGrid, 
        class LowerExtrapolationRule,  
        class UpperExtrapolationRule,  
        class... Grid1D>
class GetSplineBatchedEvaluator1D<ExecSpace, BSplines, InterpolationGrid, 
        LowerExtrapolationRule, UpperExtrapolationRule, 
        IdxRange<Grid1D...>> {
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

/// An alias describing a batched 1D spline evaluator
template<class ExecSpace, class BSplines, class InterpolationGrid, 
        class LowerExtrapolationRule,  
        class UpperExtrapolationRule,  
        class IdxRangeType>
using get_spline_batched_evaluator1d_t = typename GetSplineBatchedEvaluator1D<ExecSpace, BSplines, InterpolationGrid, 
        LowerExtrapolationRule, UpperExtrapolationRule,
        IdxRangeType>::type;


//----------------------------------------------
//  Aliases for 2D splines evaluator
//----------------------------------------------
template<class ExecSpace, 
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

template<class ExecSpace, 
        class BSplinesX1, 
        class BSplinesX2, 
        class InterpolationGridX1, 
        class InterpolationGridX2, 
        class LowerExtrapolationRuleX1,  
        class UpperExtrapolationRuleX1,  
        class LowerExtrapolationRuleX2,  
        class UpperExtrapolationRuleX2,  
        class... Grid1D>
class GetSplineBatchedEvaluator2D<ExecSpace, 
        BSplinesX1, BSplinesX2, 
        InterpolationGridX1, InterpolationGridX2, 
        LowerExtrapolationRuleX1, UpperExtrapolationRuleX1, 
        LowerExtrapolationRuleX2, UpperExtrapolationRuleX2, 
        IdxRange<Grid1D...>> {
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

/// An alias describing a batched 2D spline evaluator
template<class ExecSpace,
        class BSplinesX1, 
        class BSplinesX2, 
        class InterpolationGridX1, 
        class InterpolationGridX2, 
        class LowerExtrapolationRuleX1,  
        class UpperExtrapolationRuleX1,  
        class LowerExtrapolationRuleX2,  
        class UpperExtrapolationRuleX2,  
        class IdxRangeType>
using get_spline_batched_evaluator2d_t = typename GetSplineBatchedEvaluator2D<ExecSpace, 
        BSplinesX1, BSplinesX2, 
        InterpolationGridX1, InterpolationGridX2, 
        LowerExtrapolationRuleX1, UpperExtrapolationRuleX1, 
        LowerExtrapolationRuleX2, UpperExtrapolationRuleX2, 
        IdxRangeType>::type;

