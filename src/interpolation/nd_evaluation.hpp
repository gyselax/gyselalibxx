// SPDX-License-Identifier: MIT
#pragma once

#include "i_interpolation_evaluator.hpp"
#include "vector_field.hpp"

/**
 * @file nd_evaluation.hpp
 *
 * Free functions that lift the @c concepts::InterpolationEvaluator interface to operate
 * on VectorFields.  Each function corresponds to one of the methods defined by the
 * concept and applies it component-wise, dispatching to each scalar component via
 * @c ddcHelper::get<VectorDim>.
 */

namespace ndEval {

/**
 * @brief Evaluate an ND interpolation on each component of a VectorField.
 *
 * For each dimension @c VectorDim in @c VectorDims..., calls
 * @code
 *   evaluator(ddcHelper::get<VectorDim>(output), ddcHelper::get<VectorDim>(coeffs));
 * @endcode
 *
 * This is the VectorField counterpart of @c evaluator.operator()(eval, coeffs).
 *
 * @tparam Evaluator   A class satisfying @c concepts::InterpolationEvaluator.
 * @tparam VectorDims  The dimensions labelling the vector components (passed to
 *                     @c ddcHelper::get to extract the corresponding scalar field).
 *
 * @param[in]  evaluator  The ND interpolation evaluator.
 * @param[out] output     On output, the interpolated values for each vector component.
 * @param[in]  coeffs     The interpolation coefficients for each vector component.
 */
template <
        concepts::InterpolationEvaluator Evaluator,
        class CoordType,
        class ElementType,
        class CoeffIdxRange,
        class MemorySpace,
        class LayoutCoeff,
        class... VectorDims>
KOKKOS_FUNCTION Vector<ElementType, VectorIndexSet<VectorDims...>> evaluate(
        Evaluator const& evaluator,
        CoordType const& coord,
        VectorField<
                const ElementType,
                CoeffIdxRange,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutCoeff> coeffs)
{
    static_assert(detail::is_tagged_vector_v<CoordType>);
    (evaluator(coord, ddcHelper::get<VectorDims>(coeffs)), ...);
}

/**
 * @brief Evaluate an ND interpolation on each component of a VectorField.
 *
 * For each dimension @c VectorDim in @c VectorDims..., calls
 * @code
 *   evaluator(ddcHelper::get<VectorDim>(output), ddcHelper::get<VectorDim>(coeffs));
 * @endcode
 *
 * This is the VectorField counterpart of @c evaluator.operator()(eval, coeffs).
 *
 * @tparam Evaluator   A class satisfying @c concepts::InterpolationEvaluator.
 * @tparam VectorDims  The dimensions labelling the vector components (passed to
 *                     @c ddcHelper::get to extract the corresponding scalar field).
 *
 * @param[in]  evaluator  The ND interpolation evaluator.
 * @param[out] output     On output, the interpolated values for each vector component.
 * @param[in]  coeffs     The interpolation coefficients for each vector component.
 */
template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class EvalIdxRange,
        class CoeffIdxRange,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoeff,
        class... VectorDims>
void evaluate(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                EvalIdxRange,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        VectorField<
                const ElementType,
                CoeffIdxRange,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutCoeff> coeffs)
{
    (evaluator(ddcHelper::get<VectorDims>(output), ddcHelper::get<VectorDims>(coeffs)), ...);
}

/**
 * @brief Evaluate an ND interpolation at explicit coordinates on each component of a VectorField.
 *
 * For each dimension @c VectorDim in @c VectorDims..., calls
 * @code
 *   evaluator(ddcHelper::get<VectorDim>(output), coords, ddcHelper::get<VectorDim>(coeffs));
 * @endcode
 *
 * This is the VectorField counterpart of @c evaluator.operator()(eval, coords, coeffs).
 *
 * @tparam Evaluator   A class satisfying @c concepts::InterpolationEvaluator.
 * @tparam VectorDims  The dimensions labelling the vector components (passed to
 *                     @c ddcHelper::get to extract the corresponding scalar field).
 *
 * @param[in]  evaluator  The ND interpolation evaluator.
 * @param[out] output     On output, the interpolated values for each vector component.
 * @param[in]  coords     The evaluation coordinates.
 * @param[in]  coeffs     The interpolation coefficients for each vector component.
 */
template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class EvalIdxRange,
        class CoordType,
        class CoordIdxRange,
        class CoeffIdxRange,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoords,
        class LayoutCoeff,
        class... VectorDims>
void evaluate(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                EvalIdxRange,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        ConstField<CoordType, CoordIdxRange, MemorySpace, LayoutCoords> coords,
        VectorField<
                const ElementType,
                CoeffIdxRange,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutCoeff> coeffs)
{
    (evaluator(ddcHelper::get<VectorDims>(output), coords, ddcHelper::get<VectorDims>(coeffs)), ...);
}

/**
 * @brief Compute first-order partial derivatives of a scalar function as a VectorField.
 *
 * For each dimension @c VectorDim in @c VectorDims..., evaluates the first derivative of
 * the scalar function defined by @c coeffs in the direction @c VectorDim::Dual, storing
 * the result in the corresponding component of @c output:
 * @code
 *   evaluator.deriv(Idx<ddc::Deriv<typename VectorDim::Dual>>(1),
 *                   ddcHelper::get<VectorDim>(output), coeffs);
 * @endcode
 *
 * This is the VectorField counterpart of @c evaluator.deriv(deriv_order, eval, coeffs),
 * computing partial derivatives in all directions simultaneously.
 *
 * @tparam Evaluator   A class satisfying @c concepts::InterpolationEvaluator.
 * @tparam VectorDims  The dimensions labelling the vector components. Each @c VectorDim
 *                     must expose a @c Dual type identifying the continuous dimension
 *                     in which the derivative is taken.
 *
 * @param[in]  evaluator  The ND interpolation evaluator.
 * @param[out] output     On output, each component holds the first partial derivative
 *                        of the scalar function in the direction @c VectorDim::Dual.
 * @param[in]  coeffs     The scalar interpolation coefficients.
 */
template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class EvalIdxRange,
        class CoeffIdxRange,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoeff,
        class... VectorDims>
void deriv(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                EvalIdxRange,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        ConstField<ElementType, CoeffIdxRange, MemorySpace, LayoutCoeff> coeffs)
{
    static_assert((VectorDims::IS_COVARIANT && ...));
    (evaluator.deriv(
             Idx<ddc::Deriv<typename VectorDims::Dual>>(1),
             ddcHelper::get<VectorDims>(output),
             coeffs),
     ...);
}

/**
 * @brief Compute first-order partial derivatives at explicit coordinates as a VectorField.
 *
 * For each dimension @c VectorDim in @c VectorDims..., evaluates the first derivative of
 * the scalar function defined by @c coeffs in the direction @c VectorDim::Dual at positions
 * given by @c coords, storing the result in the corresponding component of @c output:
 * @code
 *   evaluator.deriv(Idx<ddc::Deriv<typename VectorDim::Dual>>(1),
 *                   ddcHelper::get<VectorDim>(output), coords, coeffs);
 * @endcode
 *
 * This is the VectorField counterpart of @c evaluator.deriv(deriv_order, eval, coords, coeffs).
 *
 * @tparam Evaluator   A class satisfying @c concepts::InterpolationEvaluator.
 * @tparam VectorDims  The dimensions labelling the vector components. Each @c VectorDim
 *                     must expose a @c Dual type identifying the continuous dimension
 *                     in which the derivative is taken.
 *
 * @param[in]  evaluator  The ND interpolation evaluator.
 * @param[out] output     On output, each component holds the first partial derivative
 *                        of the scalar function in the direction @c VectorDim::Dual.
 * @param[in]  coords     The evaluation coordinates.
 * @param[in]  coeffs     The scalar interpolation coefficients.
 */
template <
        concepts::InterpolationEvaluator Evaluator,
        class ElementType,
        class EvalIdxRange,
        class CoordType,
        class CoordIdxRange,
        class CoeffIdxRange,
        class MemorySpace,
        class LayoutOut,
        class LayoutCoords,
        class LayoutCoeff,
        class... VectorDims>
void deriv(
        Evaluator const& evaluator,
        VectorField<
                ElementType,
                EvalIdxRange,
                VectorIndexSet<VectorDims...>,
                MemorySpace,
                LayoutOut> output,
        ConstField<CoordType, CoordIdxRange, MemorySpace, LayoutCoords> coords,
        ConstField<ElementType, CoeffIdxRange, MemorySpace, LayoutCoeff> coeffs)
{
    static_assert((VectorDims::IS_COVARIANT && ...));
    (evaluator.deriv(
             Idx<ddc::Deriv<typename VectorDims::Dual>>(1),
             ddcHelper::get<VectorDims>(output),
             coords,
             coeffs),
     ...);
}

} // namespace ndEval
