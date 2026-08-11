

# File vector\_field\_evaluation.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**vector\_field\_evaluation.hpp**](vector__field__evaluation_8hpp.md)

[Go to the documentation of this file](vector__field__evaluation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "i_interpolation_evaluator.hpp"
#include "vector_field.hpp"

namespace ndEval {

template <
        concepts::InterpolationEvaluator Evaluator,
        class CoordType,
        class ElementType,
        class IdxRangeCoeff,
        class MemorySpace,
        class LayoutCoeff,
        class... VectorDims>
KOKKOS_FUNCTION Vector<ElementType, VectorDims...> evaluate(
        Evaluator const& evaluator,
        CoordType const& coord,
        VectorField<
                const ElementType,
                IdxRangeCoeff,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutCoeff> coeffs)
{
    static_assert(ddc::detail::is_tagged_vector_v<CoordType>);
    return Vector<ElementType, VectorDims...>(
            evaluator(coord, ddcHelper::get<VectorDims>(coeffs))...);
}

template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class IdxRangeEval,
        class IdxRangeCoeff,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoeff,
        class... VectorDims>
void evaluate(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                IdxRangeEval,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        VectorField<
                const ElementType,
                IdxRangeCoeff,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutCoeff> coeffs)
{
    (evaluator(ddcHelper::get<VectorDims>(output), ddcHelper::get<VectorDims>(coeffs)), ...);
}

template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class IdxRangeEval,
        class CoordType,
        class IdxRangeCoord,
        class IdxRangeCoeff,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoords,
        class LayoutCoeff,
        class... VectorDims>
void evaluate(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                IdxRangeEval,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        ConstField<CoordType, IdxRangeCoord, MemorySpace, LayoutCoords> coords,
        VectorField<
                const ElementType,
                IdxRangeCoeff,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutCoeff> coeffs)
{
    (evaluator(ddcHelper::get<VectorDims>(output), coords, ddcHelper::get<VectorDims>(coeffs)),
     ...);
}

template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class IdxRangeEval,
        class IdxRangeCoeff,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoeff,
        class... VectorDims>
void deriv(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                IdxRangeEval,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        ConstField<ElementType, IdxRangeCoeff, MemorySpace, LayoutCoeff> coeffs)
{
    static_assert((VectorDims::IS_COVARIANT && ...));
    (evaluator
             .deriv(Idx<ddc::Deriv<typename VectorDims::Dual>>(1),
                    ddcHelper::get<VectorDims>(output),
                    coeffs),
     ...);
}

template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class IdxRangeEval,
        class CoordType,
        class IdxRangeCoord,
        class IdxRangeCoeff,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoords,
        class LayoutCoeff,
        class... VectorDims>
void deriv(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                IdxRangeEval,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        ConstField<CoordType, IdxRangeCoord, MemorySpace, LayoutCoords> coords,
        ConstField<ElementType, IdxRangeCoeff, MemorySpace, LayoutCoeff> coeffs)
{
    static_assert((VectorDims::IS_COVARIANT && ...));
    (evaluator
             .deriv(Idx<ddc::Deriv<typename VectorDims::Dual>>(1),
                    ddcHelper::get<VectorDims>(output),
                    coords,
                    coeffs),
     ...);
}

} // namespace ndEval
```


