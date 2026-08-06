

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
#include "nd_identity_interpolation_builder.hpp"
#include "nd_lagrange_evaluator.hpp"

namespace detail {

template <
        class ExecSpace,
        class Basis,
        class InterpGrid,
        class ExtrapRules,
        class DataType = double>
class LagrangeInterpolator
{
    using continuous_dimension_type = typename InterpGrid::continuous_dimension_type;

    using MinExtrapolationRule = ddc::type_seq_element_t<0, ExtrapRules>;
    using MaxExtrapolationRule = ddc::type_seq_element_t<1, ExtrapRules>;

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
                    MinExtrapolationRule,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type>>);
    static_assert(
            is_periodic
            == std::is_same_v<
                    MaxExtrapolationRule,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type>>);

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
            (is_extrapolation_rule_auto_constructible_v<MinExtrapolationRule, CoeffGridType, DataType, Basis>)&&(
                    is_extrapolation_rule_auto_constructible_v<
                            MaxExtrapolationRule,
                            CoeffGridType,
                            DataType,
                            Basis>))
        : m_min_extrapolation(
                get_extrapolation<MinExtrapolationRule, CoeffGridType, DataType, Basis>(
                        Extremity::FRONT))
        , m_max_extrapolation(
                  get_extrapolation<MaxExtrapolationRule, CoeffGridType, DataType, Basis>(
                          Extremity::BACK))
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

template <
        class ExecSpace,
        class IdxRangeBasis,
        class IdxRangeInterpGrid,
        class ExtrapRulesSeq,
        class DataType = double>
class NDLagrangeInterpolator;

template <class ExecSpace, class... Basis, class... Grid1D, class... ExtrapRules, class DataType>
class NDLagrangeInterpolator<
        ExecSpace,
        IdxRange<Basis...>,
        IdxRange<Grid1D...>,
        ddc::detail::TypeSeq<ExtrapRules...>,
        DataType>
{
    static_assert(sizeof...(Basis) == sizeof...(Grid1D));
    static_assert(sizeof...(Basis) == sizeof...(ExtrapRules));
    static_assert(sizeof...(Basis) > 0);
    static_assert((is_lagrange_basis_v<Basis> && ...));

public:
    using memory_space = typename ExecSpace::memory_space;

private:
    template <class Rule>
    using MinRule = ddc::type_seq_element_t<0, Rule>;

    template <class Rule>
    using MaxRule = ddc::type_seq_element_t<1, Rule>;

public:
    using BuilderType = NDIdentityInterpolationBuilder<
            ExecSpace,
            memory_space,
            DataType,
            IdxRange<Grid1D...>,
            IdxRange<Basis...>>;

    using EvaluatorType = NDLagrangeEvaluator<LagrangeEvaluator<
            ExecSpace,
            memory_space,
            DataType,
            Basis,
            Grid1D,
            MinRule<ExtrapRules>,
            MaxRule<ExtrapRules>>...>;

private:
    template <class B>
    using CoeffGrid = ddc::type_seq_element_t<
            ddc::type_seq_rank_v<B, ddc::detail::TypeSeq<Basis...>>,
            ddc::to_type_seq_t<typename BuilderType::coeff_idx_range_type>>;

    static_assert(
            ((Basis::is_periodic()
              == std::is_same_v<
                      MinRule<ExtrapRules>,
                      ddc::PeriodicExtrapolationRule<
                              typename Basis::continuous_dimension_type>>)&&...),
            "PeriodicExtrapolationRule has to be used if and only if the dimension is "
            "periodic");
    static_assert(
            ((Basis::is_periodic()
              == std::is_same_v<
                      MaxRule<ExtrapRules>,
                      ddc::PeriodicExtrapolationRule<
                              typename Basis::continuous_dimension_type>>)&&...),
            "PeriodicExtrapolationRule has to be used if and only if the dimension is "
            "periodic");

public:
    static constexpr std::size_t rank()
    {
        return sizeof...(Grid1D);
    }

private:
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    explicit NDLagrangeInterpolator(IdxRange<Grid1D...> idx_range = IdxRange<Grid1D...> {}) requires(
            (is_extrapolation_rule_auto_constructible_v<
                     MinRule<ExtrapRules>,
                     CoeffGrid<Basis>,
                     DataType,
                     Basis> && ...)
            && (is_extrapolation_rule_auto_constructible_v<
                        MaxRule<ExtrapRules>,
                        CoeffGrid<Basis>,
                        DataType,
                        Basis> && ...))
        : m_evaluator(LagrangeEvaluator<
                      ExecSpace,
                      memory_space,
                      DataType,
                      Basis,
                      Grid1D,
                      MinRule<ExtrapRules>,
                      MaxRule<ExtrapRules>>(
                get_extrapolation<MinRule<ExtrapRules>, CoeffGrid<Basis>, DataType, Basis>(
                        Extremity::FRONT),
                get_extrapolation<MaxRule<ExtrapRules>, CoeffGrid<Basis>, DataType, Basis>(
                        Extremity::BACK))...)
    {
    }

    explicit NDLagrangeInterpolator(
            IdxRange<Grid1D...> idx_range,
            std::pair<MinRule<ExtrapRules>, MaxRule<ExtrapRules>> const&... extrapolation_rules)
        : m_evaluator(LagrangeEvaluator<
                      ExecSpace,
                      memory_space,
                      DataType,
                      Basis,
                      Grid1D,
                      MinRule<ExtrapRules>,
                      MaxRule<ExtrapRules>>(
                extrapolation_rules.first,
                extrapolation_rules.second)...)
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

template <
        class ExecSpace,
        class DataType,
        class IdxRangeBasis,
        class IdxRangeInterpGrid,
        class... ExtrapRulesSeq>
struct LagrangeInterpolatorResolver;

template <class ExecSpace, class DataType, class Basis, class Grid1D, class ExtrapRules>
struct LagrangeInterpolatorResolver<
        ExecSpace,
        DataType,
        IdxRange<Basis>,
        IdxRange<Grid1D>,
        ExtrapRules>
{
    using type = detail::LagrangeInterpolator<
            ExecSpace,
            Basis,
            Grid1D,
            extrapolation_rule_t<
                    ExtrapRules,
                    typename IdentityInterpolationBuilder<
                            ExecSpace,
                            typename ExecSpace::memory_space,
                            DataType,
                            Grid1D,
                            Basis>::basis_domain_type,
                    DataType,
                    Basis>,
            DataType>;
};

template <
        class ExecSpace,
        class DataType,
        class BasisHead,
        class... Basis,
        class Grid1DHead,
        class... Grid1D,
        class ExtrapRulesHead,
        class... ExtrapRules>
struct LagrangeInterpolatorResolver<
        ExecSpace,
        DataType,
        IdxRange<BasisHead, Basis...>,
        IdxRange<Grid1DHead, Grid1D...>,
        ExtrapRulesHead,
        ExtrapRules...>
{
    using type = detail::NDLagrangeInterpolator<
            ExecSpace,
            IdxRange<BasisHead, Basis...>,
            IdxRange<Grid1DHead, Grid1D...>,
            ddc::detail::TypeSeq<
                    extrapolation_rule_t<
                            ExtrapRulesHead,
                            typename BasisHead::template Impl<
                                    BasisHead,
                                    typename ExecSpace::memory_space>::knot_grid,
                            DataType,
                            BasisHead>,
                    extrapolation_rule_t<
                            ExtrapRules,
                            typename Basis::template Impl<Basis, typename ExecSpace::memory_space>::
                                    knot_grid,
                            DataType,
                            Basis>...>,
            DataType>;
};
} // namespace detail

template <
        class ExecSpace,
        class DataType,
        class IdxRangeBasis,
        class IdxRangeInterpGrid,
        class... ExtrapRules>
using LagrangeInterpolator = typename detail::LagrangeInterpolatorResolver<
        ExecSpace,
        DataType,
        IdxRangeBasis,
        IdxRangeInterpGrid,
        ExtrapRules...>::type;
```


