

# File i\_interpolation\_evaluator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**i\_interpolation\_evaluator.hpp**](i__interpolation__evaluator_8hpp.md)

[Go to the documentation of this file](i__interpolation__evaluator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <type_traits>

#include "ddc_aliases.hpp"

template <class Evaluator>
struct InterpolationEvaluatorTraits
{
    using data_type = typename Evaluator::data_type;

    using evaluation_idx_range_type = typename Evaluator::evaluation_idx_range_type;

    using coeff_grid_type = typename Evaluator::coeff_grid_type;

    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_idx_range_type<
                    BatchedInterpolationIdxRange>;

    template <class D>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_coeff_idx_range_type<D>;
};

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
    using data_type = double;

    using evaluation_idx_range_type = typename Evaluator::evaluation_domain_type;

    using coeff_grid_type = typename Evaluator::bsplines_type;

    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_domain_type<
                    BatchedInterpolationIdxRange>;

    template <class BatchedInterpolationIdxRange>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_spline_domain_type<BatchedInterpolationIdxRange>;
};

namespace concepts {

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
```


