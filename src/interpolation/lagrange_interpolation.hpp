// SPDX-License-Identifier: MIT
#pragma once

#include "extrapolation_rule_choice.hpp"
#include "identity_interpolation_builder.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "lagrange_evaluator.hpp"

template <
        class ExecSpace,
        class Basis,
        class InterpGrid,
        ExtrapolationRule MinExtrapRule,
        ExtrapolationRule MaxExtrapRule,
        ddc::BoundCond MinBound = ddc::BoundCond::GREVILLE,
        ddc::BoundCond MaxBound = ddc::BoundCond::GREVILLE,
        class DataType = double>
class LagrangeInterpolator
{
    using continuous_dimension_type = typename InterpGrid::continuous_dimension_type;

    static constexpr bool is_periodic = continuous_dimension_type::PERIODIC;

    static_assert(is_periodic == (MinBound == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic == (MaxBound == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic == (MinExtrapRule == ExtrapolationRule::PERIODIC));
    static_assert(is_periodic == (MaxExtrapRule == ExtrapolationRule::PERIODIC));

    static_assert(is_lagrange_basis_v<Basis>);

public:
    using BuilderType = IdentityInterpolationBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
            InterpGrid,
            Basis>;

    using CoeffGridType = typename BuilderType::basis_domain_type;

    using EvaluatorType = LagrangeEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
            Basis,
            InterpGrid,
            extrapolation_rule_t<MinExtrapRule, CoeffGridType>,
            extrapolation_rule_t<MaxExtrapRule, CoeffGridType>>;

private:
    extrapolation_rule_t<MinExtrapRule, CoeffGridType> m_min_extrapolation;
    extrapolation_rule_t<MaxExtrapRule, CoeffGridType> m_max_extrapolation;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    LagrangeInterpolator()
        : m_min_extrapolation(
                get_extrapolation<MinExtrapRule, CoeffGridType, Basis>(Extremity::FRONT))
        , m_max_extrapolation(
                  get_extrapolation<MaxExtrapRule, CoeffGridType, Basis>(Extremity::BACK))
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
