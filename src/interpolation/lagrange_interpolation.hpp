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
        ddc::BoundCond MinBound,
        ddc::BoundCond MaxBound,
        ExtrapolationRule MinExtrapRule,
        ExtrapolationRule MaxExtrapRule,
        class DataType = double>
class LagrangeInterpolator
{
    // Constant not yet implemented for Lagrange
    static_assert(MinExtrapRule != ExtrapolationRule::CONSTANT);
    static_assert(MaxExtrapRule != ExtrapolationRule::CONSTANT);

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

    using EvaluatorType = LagrangeEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
            Basis,
            InterpGrid,
            details::extrapolation_rule_t<MinExtrapRule, Basis>,
            details::extrapolation_rule_t<MaxExtrapRule, Basis>>;

private:
    details::extrapolation_rule_t<MinExtrapRule, Basis> m_min_extrapolation;
    details::extrapolation_rule_t<MaxExtrapRule, Basis> m_max_extrapolation;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    LagrangeInterpolator(IdxRange<InterpGrid> idx_range)
        : m_min_extrapolation(details::get_extrapolation<MinExtrapRule, Basis>(Extremity::FRONT))
        , m_max_extrapolation(details::get_extrapolation<MaxExtrapRule, Basis>(Extremity::BACK))
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
