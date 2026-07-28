

# File lagrange\_interpolation.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**lagrange\_interpolation.hpp**](lagrange__interpolation_8hpp.md)

[Go to the documentation of this file](lagrange__interpolation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <utility>

#include "extrapolation_rule_choice.hpp"
#include "identity_interpolation_builder.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "lagrange_evaluator.hpp"

template <
        class ExecSpace,
        class Basis,
        class InterpGrid,
        class MinExtrapRule,
        class MaxExtrapRule,
        class DataType = double>
class LagrangeInterpolator
{
    using continuous_dimension_type = typename InterpGrid::continuous_dimension_type;

    static constexpr bool is_periodic = continuous_dimension_type::PERIODIC;

    static_assert(is_lagrange_basis_v<Basis>);

public:
    using BuilderType = IdentityInterpolationBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
            InterpGrid,
            Basis>;

    using CoeffGridType = typename BuilderType::basis_domain_type;

private:
    static_assert(
            is_periodic
            == std::is_same_v<
                    extrapolation_rule_t<MinExtrapRule, CoeffGridType, DataType>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type>>);
    static_assert(
            is_periodic
            == std::is_same_v<
                    extrapolation_rule_t<MaxExtrapRule, CoeffGridType, DataType>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type>>);

    using MinExtrapolationRule = extrapolation_rule_t<MinExtrapRule, CoeffGridType, DataType>;
    using MaxExtrapolationRule = extrapolation_rule_t<MaxExtrapRule, CoeffGridType, DataType>;

public:
    using EvaluatorType = LagrangeEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
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
    explicit LagrangeInterpolator(IdxRange<InterpGrid> idx_range = IdxRange<InterpGrid> {}) requires(
            is_extrapolation_rule_auto_constructible_v<MinExtrapRule, CoeffGridType, DataType>&&
                    is_extrapolation_rule_auto_constructible_v<
                            MaxExtrapRule,
                            CoeffGridType,
                            DataType>)
        : m_min_extrapolation(
                get_extrapolation<MinExtrapRule, CoeffGridType, DataType, Basis>(Extremity::FRONT))
        , m_max_extrapolation(
                  get_extrapolation<MaxExtrapRule, CoeffGridType, DataType, Basis>(Extremity::BACK))
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

    explicit LagrangeInterpolator(
            IdxRange<InterpGrid> idx_range,
            MinExtrapolationRule min_extrapolation_rule,
            MaxExtrapolationRule max_extrapolation_rule)
        : m_min_extrapolation(std::move(min_extrapolation_rule))
        , m_max_extrapolation(std::move(max_extrapolation_rule))
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


