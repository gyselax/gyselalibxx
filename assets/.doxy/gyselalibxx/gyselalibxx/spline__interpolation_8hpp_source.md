

# File spline\_interpolation.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**spline\_interpolation.hpp**](spline__interpolation_8hpp.md)

[Go to the documentation of this file](spline__interpolation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <utility>

#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "extrapolation_rule_choice.hpp"

template <
        class ExecSpace,
        class Basis,
        class InterpGrid,
        class MinExtrapRule,
        class MaxExtrapRule,
        ddc::SplineBuilderClosure MinBound,
        ddc::SplineBuilderClosure MaxBound,
        ddc::SplineSolver Solver = ddc::SplineSolver::LAPACK>
class SplineInterpolator
{
private:
    using continuous_dimension_type = typename InterpGrid::continuous_dimension_type;

    static constexpr bool is_periodic = continuous_dimension_type::PERIODIC;

    static_assert(is_periodic == (MinBound == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic == (MaxBound == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic
            == std::is_same_v<
                    extrapolation_rule_t<MinExtrapRule, Basis, double>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type>>);
    static_assert(
            is_periodic
            == std::is_same_v<
                    extrapolation_rule_t<MaxExtrapRule, Basis, double>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type>>);

    static_assert(is_spline_basis_v<Basis>);

    using MinExtrapolationRule = extrapolation_rule_t<MinExtrapRule, Basis, double>;
    using MaxExtrapolationRule = extrapolation_rule_t<MaxExtrapRule, Basis, double>;

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
            MinExtrapolationRule,
            MaxExtrapolationRule>;

    static constexpr std::size_t rank()
    {
        return 1;
    }

private:
    MinExtrapolationRule m_min_extrapolation;
    MaxExtrapolationRule m_max_extrapolation;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    explicit SplineInterpolator(std::string const& label, IdxRange<InterpGrid> idx_range) requires(
            is_extrapolation_rule_auto_constructible_v<MinExtrapRule, Basis, double>&&
                    is_extrapolation_rule_auto_constructible_v<MaxExtrapRule, Basis, double>)
        : m_min_extrapolation(get_extrapolation<MinExtrapRule, Basis, double>(Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<MaxExtrapRule, Basis, double>(Extremity::BACK))
        , m_builder(label, idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }
    explicit SplineInterpolator(IdxRange<InterpGrid> idx_range) requires(
            is_extrapolation_rule_auto_constructible_v<MinExtrapRule, Basis, double>&&
                    is_extrapolation_rule_auto_constructible_v<MaxExtrapRule, Basis, double>)
        : m_min_extrapolation(get_extrapolation<MinExtrapRule, Basis, double>(Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<MaxExtrapRule, Basis, double>(Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

    explicit SplineInterpolator(
            std::string const& label,
            IdxRange<InterpGrid> idx_range,
            MinExtrapolationRule min_extrapolation_rule,
            MaxExtrapolationRule max_extrapolation_rule)
        : m_min_extrapolation(std::move(min_extrapolation_rule))
        , m_max_extrapolation(std::move(max_extrapolation_rule))
        , m_builder(label, idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

    explicit SplineInterpolator(
            IdxRange<InterpGrid> idx_range,
            MinExtrapolationRule min_extrapolation_rule,
            MaxExtrapolationRule max_extrapolation_rule)
        : m_min_extrapolation(std::move(min_extrapolation_rule))
        , m_max_extrapolation(std::move(max_extrapolation_rule))
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


