

# File spline\_interpolation.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**spline\_interpolation.hpp**](spline__interpolation_8hpp.md)

[Go to the documentation of this file](spline__interpolation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "extrapolation_rule_choice.hpp"

template <
        class ExecSpace,
        class Basis,
        class InterpGrid,
        ExtrapolationRule MinExtrapRule,
        ExtrapolationRule MaxExtrapRule,
        ddc::BoundCond MinBound,
        ddc::BoundCond MaxBound,
        ddc::SplineSolver Solver = ddc::SplineSolver::LAPACK>
class SplineInterpolator
{
private:
    using continuous_dimension_type = typename InterpGrid::continuous_dimension_type;

    static constexpr bool is_periodic = continuous_dimension_type::PERIODIC;

    static_assert(is_periodic == (MinBound == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic == (MaxBound == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic == (MinExtrapRule == ExtrapolationRule::PERIODIC));
    static_assert(is_periodic == (MaxExtrapRule == ExtrapolationRule::PERIODIC));

    static_assert(is_spline_basis_v<Basis>);

public:
    using BuilderType = ddc::SplineBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis,
            InterpGrid,
            MinBound,
            MaxBound,
            Solver>;

    using EvaluatorType = ddc::SplineEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis,
            InterpGrid,
            extrapolation_rule_t<MinExtrapRule, Basis>,
            extrapolation_rule_t<MaxExtrapRule, Basis>>;

    static constexpr std::size_t rank()
    {
        return 1;
    }

private:
    extrapolation_rule_t<MinExtrapRule, Basis> m_min_extrapolation;
    extrapolation_rule_t<MaxExtrapRule, Basis> m_max_extrapolation;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    explicit SplineInterpolator(IdxRange<InterpGrid> idx_range)
        : m_min_extrapolation(get_extrapolation<MinExtrapRule, Basis>(Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<MaxExtrapRule, Basis>(Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    EvaluatorType const& get_evaluator() const
    {
        return m_evaluator;
    }
};
```


