// SPDX-License-Identifier: MIT
#pragma once

#include <type_traits>

#include "ddc_aliases.hpp"

/**
 * @brief A traits struct for accessing type aliases of an interpolation evaluator.
 *
 * The primary template delegates to the evaluator's own type aliases. It handles
 * the common case where an evaluator already uses the convention names (e.g.
 * LagrangeEvaluator exposes @c evaluation_idx_range_type and
 * @c batched_coeff_idx_range_type directly).
 *
 * Specialise this struct to adapt external evaluators whose alias names differ
 * (e.g. ddc::SplineEvaluator).
 *
 * @tparam Evaluator The interpolation evaluator type.
 */
template <class Evaluator>
struct InterpolationEvaluatorTraits
{
    /// @brief The data type that the data is saved on.
    using data_type = typename Evaluator::data_type;

    /// @brief The 1D index range for the evaluation mesh.
    using evaluation_idx_range_type = typename Evaluator::evaluation_idx_range_type;

    /// @brief The discrete dimension for the interpolation coefficients.
    using coeff_grid_type = typename Evaluator::coeff_grid_type;

    /// @brief Batched index range for the evaluation
    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_idx_range_type<
                    BatchedInterpolationIdxRange>;

    /// @brief Batched domain with the evaluation grid replaced by coeff_grid_type.
    template <class D>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_coeff_idx_range_type<D>;
};

/**
 * @brief Specialisation of InterpolationEvaluatorTraits for ddc::SplineEvaluator.
 *
 * ddc::SplineEvaluator uses different alias names from the InterpolationEvaluator
 * convention. This specialisation provides the mapping so that ddc::SplineEvaluator
 * can be used directly as an InterpolationEvaluator.
 *
 * Mapping:
 *   evaluation_discrete_dimension_type -> (defines evaluation_idx_range_type)
 *   evaluation_domain_type             -> evaluation_idx_range_type
 *   bsplines_type                      -> coeff_grid_type
 *   batched_spline_domain_type<D>      -> batched_coeff_idx_range_type<D>
 */
template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class EvaluationDDim,
        class LowerExtrapolationRule,
        class UpperExtrapolationRule>
struct InterpolationEvaluatorTraits<ddc::SplineEvaluator<
        ExecSpace,
        MemorySpace,
        BSplines,
        EvaluationDDim,
        LowerExtrapolationRule,
        UpperExtrapolationRule>>
{
private:
    using Evaluator = ddc::SplineEvaluator<
            ExecSpace,
            MemorySpace,
            BSplines,
            EvaluationDDim,
            LowerExtrapolationRule,
            UpperExtrapolationRule>;

public:
    /// @brief The data type that the data is saved on.
    using data_type = double;

    /// @brief The 1D index range for the evaluation mesh.
    using evaluation_idx_range_type = typename Evaluator::evaluation_domain_type;

    /// @brief The discrete dimension for the B-spline coefficients.
    using coeff_grid_type = typename Evaluator::bsplines_type;

    /// @brief Batched index range for the evaluation
    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_domain_type<
                    BatchedInterpolationIdxRange>;

    /// @brief Batched domain with the evaluation grid replaced by BSplines.
    template <class BatchedInterpolationIdxRange>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_spline_domain_type<BatchedInterpolationIdxRange>;
};

namespace concepts {

/**
 * @brief A concept describing an interpolation evaluator.
 *
 * An interpolation evaluator is a callable that takes a field of interpolation
 * coefficients (on a basis domain) and produces function values on an evaluation
 * mesh, optionally at user-supplied coordinates.
 *
 * Type information is accessed through InterpolationEvaluatorTraits<Evaluator>, which
 * has a primary template delegating to the evaluator's own aliases and can be
 * specialised for external evaluators (e.g. ddc::SplineEvaluator).
 *
 * The template alias and method requirements are verified using evaluation_idx_range_type
 * as a representative domain.
 */
template <class Evaluator>
concept InterpolationEvaluator = requires
{
    typename Evaluator::exec_space;
    typename Evaluator::memory_space;
    typename Evaluator::continuous_dimension_type;
    typename InterpolationEvaluatorTraits<Evaluator>::data_type;
    typename InterpolationEvaluatorTraits<Evaluator>::evaluation_idx_range_type;
    typename InterpolationEvaluatorTraits<Evaluator>::coeff_grid_type;
    typename Evaluator::lower_extrapolation_rule_type;
    typename Evaluator::upper_extrapolation_rule_type;
    // Verify the template alias can be instantiated with a concrete domain
    typename InterpolationEvaluatorTraits<Evaluator>::template batched_coeff_idx_range_type<
            typename InterpolationEvaluatorTraits<Evaluator>::evaluation_idx_range_type>;
}
&&requires(
        Evaluator const& e,
        Field<typename InterpolationEvaluatorTraits<Evaluator>::data_type,
              typename InterpolationEvaluatorTraits<Evaluator>::evaluation_idx_range_type,
              typename Evaluator::memory_space> eval,
        ConstField<
                typename InterpolationEvaluatorTraits<Evaluator>::data_type,
                typename InterpolationEvaluatorTraits<Evaluator>::
                        template batched_coeff_idx_range_type<typename InterpolationEvaluatorTraits<
                                Evaluator>::evaluation_idx_range_type>,
                typename Evaluator::memory_space> coeffs,
        ConstField<
                Coord<typename Evaluator::continuous_dimension_type>,
                typename InterpolationEvaluatorTraits<Evaluator>::evaluation_idx_range_type,
                typename Evaluator::memory_space> coords)
{
    {e(eval, coeffs)};
    {e(eval, coords, coeffs)};
}
&&requires(
        Evaluator const& e,
        ConstField<
                typename InterpolationEvaluatorTraits<Evaluator>::data_type,
                typename InterpolationEvaluatorTraits<Evaluator>::
                        template batched_coeff_idx_range_type<typename InterpolationEvaluatorTraits<
                                Evaluator>::evaluation_idx_range_type>,
                typename Evaluator::memory_space> coeffs,
        Coord<typename Evaluator::continuous_dimension_type> coord)
{
    {
        e(coord, coeffs)
        } -> std::same_as<typename InterpolationEvaluatorTraits<Evaluator>::data_type>;
};

} // namespace concepts
