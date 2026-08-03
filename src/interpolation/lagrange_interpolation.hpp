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

/**
 * @brief An owning interpolation object that bundles a Lagrange builder and evaluator.
 *
 * LagrangeInterpolator constructs and owns a matching IdentityInterpolationBuilder and
 * LagrangeEvaluator for a given dimension. It satisfies the concepts::Interpolation
 * concept and is the recommended way to create a Lagrange interpolation for use with
 * advection operators and similar algorithms.
 *
 * The builder is an identity operation: it passes function values on the interpolation
 * grid directly as coefficients to the evaluator, which then performs local polynomial
 * reconstruction via the Lagrange basis.
 *
 * @tparam ExecSpace     The Kokkos execution space used for computations.
 * @tparam Basis         The Lagrange basis type (uniform or non-uniform).
 * @tparam InterpGrid    The discrete grid on which function values are provided.
 * @tparam ExtrapRules   A ddc::detail::TypeSeq<MinExtrapolationRule, MaxExtrapolationRule> pairing the
 *                       extrapolation rules applied below/above the boundary. Where
 *                       MinExtrapolationRule and MaxExtrapolationRule are extrapolation rule classes.
 * @tparam DataType      The floating-point type of the function values (default: double).
 */
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
    /// @brief The IdentityInterpolationBuilder type built from the template parameters.
    using BuilderType = IdentityInterpolationBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
            InterpGrid,
            Basis>;

    /// @brief The discrete grid type used for the Lagrange coefficients (the Lagrange basis grid).
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
    /// @brief The LagrangeEvaluator type built from the template parameters.
    using EvaluatorType = LagrangeEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
            Basis,
            InterpGrid,
            MinExtrapolationRule,
            MaxExtrapolationRule>;

    /// @brief The number of interpolation dimensions.
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
    /**
     * @brief Construct a LagrangeInterpolator.
     *
     * The extrapolation rules are initialised from the discrete space of @c Basis,
     * so the corresponding ddc discrete space must be initialised before construction.
     * No index range is required because the identity builder needs none.
     * This overload is only available when both extrapolation rules can be built
     * automatically (they are default-constructible, or the tag ExtrapolationRule::Constant
     * is used) - otherwise use the overload that takes the extrapolation rules explicitly.
     *
     * @param idx_range The index range on which the interpolator will act. This is
     *                  unused but is included to match the SplineInterpolator interface.
     */
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

    /**
     * @brief Construct a LagrangeInterpolator, specifying the extrapolation rules explicitly.
     *
     * Use this overload when the chosen extrapolation rule cannot be built automatically,
     * e.g. a custom extrapolation rule that is not default-constructible and is not
     * ExtrapolationRule::Constant.
     *
     * @param idx_range The index range on which the interpolator will act. This is
     *                  unused but is included to match the SplineInterpolator interface.
     * @param min_extrapolation_rule The extrapolation rule to use below the lower boundary.
     * @param max_extrapolation_rule The extrapolation rule to use above the upper boundary.
     */
    explicit LagrangeInterpolator(
            IdxRange<InterpGrid> idx_range,
            MinExtrapolationRule min_extrapolation_rule,
            MaxExtrapolationRule max_extrapolation_rule)
        : m_min_extrapolation(std::move(min_extrapolation_rule))
        , m_max_extrapolation(std::move(max_extrapolation_rule))
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

    /**
     * @brief Return a const reference to the owned identity builder.
     * @return The BuilderType instance.
     */
    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    /**
     * @brief Return a const reference to the owned Lagrange evaluator.
     * @return The EvaluatorType instance.
     */
    EvaluatorType const& get_evaluator() const
    {
        return m_evaluator;
    }
};

/**
 * @brief An owning interpolation object that bundles an ND Lagrange builder and evaluator.
 *
 * NDLagrangeInterpolator constructs and owns a matching NDIdentityInterpolationBuilder and
 * NDLagrangeEvaluator for a tensor-product grid of dimensions. It satisfies the
 * concepts::Interpolation concept and is the ND generalisation of LagrangeInterpolator.
 *
 * The builder is an identity operation: it passes function values on the interpolation
 * mesh directly as coefficients to the evaluator, which then performs local polynomial
 * reconstruction via a tensor product of 1D Lagrange bases.
 *
 * @tparam ExecSpace     The Kokkos execution space used for computations.
 * @tparam IdxRangeBasis The ND index range for the Lagrange basis types, of the form
 *                       IdxRange<Basis1, ..., BasisN>, one per interpolation dimension.
 * @tparam IdxRangeInterpGrid The ND index range for the interpolation mesh, of the form
 *                       IdxRange<Grid1, ..., GridN>, in the same order as IdxRangeBasis.
 * @tparam ExtrapRulesSeq A ddc::detail::TypeSeq<ExtrapRules1, ..., ExtrapRulesN> pairing,
 *                       for each dimension, the extrapolation rules applied below/above
 *                       the boundary (each ExtrapRulesI is itself a
 *                       ddc::detail::TypeSeq<MinExtrapolationRuleI, MaxExtrapolationRuleI>).
 * @tparam DataType      The floating-point type of the function values (default: double).
 */
template <
        class ExecSpace,
        class IdxRangeBasis,
        class IdxRangeInterpGrid,
        class ExtrapRulesSeq,
        class DataType = double>
class NDLagrangeInterpolator;

/// The implementation of NDLagrangeInterpolator. This is separate to allow variadic packs.
template <class ExecSpace, class... Basis, class... Grid, class... ExtrapRules, class DataType>
class NDLagrangeInterpolator<
        ExecSpace,
        IdxRange<Basis...>,
        IdxRange<Grid...>,
        ddc::detail::TypeSeq<ExtrapRules...>,
        DataType>
{
    static_assert(sizeof...(Basis) == sizeof...(Grid));
    static_assert(sizeof...(Basis) == sizeof...(ExtrapRules));
    static_assert(sizeof...(Basis) > 0);
    static_assert((is_lagrange_basis_v<Basis> && ...));

public:
    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = typename ExecSpace::memory_space;

private:
    template <class Rule>
    using MinOf = ddc::type_seq_element_t<0, Rule>;

    template <class Rule>
    using MaxOf = ddc::type_seq_element_t<1, Rule>;

    template <class B>
    using CoeffGridOf = typename B::template Impl<B, memory_space>::knot_grid;

    static_assert(
            ((Basis::is_periodic()
              == std::is_same_v<
                      MinOf<ExtrapRules>,
                      ddc::PeriodicExtrapolationRule<
                              typename Basis::continuous_dimension_type>>)&&...),
            "PeriodicExtrapolationRule has to be used if and only if the dimension is "
            "periodic");
    static_assert(
            ((Basis::is_periodic()
              == std::is_same_v<
                      MaxOf<ExtrapRules>,
                      ddc::PeriodicExtrapolationRule<
                              typename Basis::continuous_dimension_type>>)&&...),
            "PeriodicExtrapolationRule has to be used if and only if the dimension is "
            "periodic");

public:
    /// @brief The NDIdentityInterpolationBuilder type built from the template parameters.
    using BuilderType = NDIdentityInterpolationBuilder<
            ExecSpace,
            memory_space,
            DataType,
            IdxRange<Grid...>,
            IdxRange<Basis...>>;

    /// @brief The NDLagrangeEvaluator type built from the template parameters.
    using EvaluatorType = NDLagrangeEvaluator<LagrangeEvaluator<
            ExecSpace,
            memory_space,
            DataType,
            Basis,
            Grid,
            MinOf<ExtrapRules>,
            MaxOf<ExtrapRules>>...>;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return sizeof...(Grid);
    }

private:
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    /**
     * @brief Construct an NDLagrangeInterpolator.
     *
     * The extrapolation rules are initialised from the discrete spaces of the Lagrange
     * bases, so the corresponding ddc discrete spaces must be initialised before
     * construction. No index range is required because the identity builder needs none.
     * This overload is only available when every dimension's extrapolation rules can be
     * built automatically (see is_extrapolation_rule_auto_constructible_v) - otherwise use
     * the overload that takes the extrapolation rules explicitly.
     *
     * @param idx_range The index range on which the interpolator will act. This is
     *                  unused but is included to match the SplineInterpolator interface.
     */
    explicit NDLagrangeInterpolator(IdxRange<Grid...> idx_range = IdxRange<Grid...> {}) requires(
            (is_extrapolation_rule_auto_constructible_v<
                     MinOf<ExtrapRules>,
                     CoeffGridOf<Basis>,
                     DataType,
                     Basis> && ...)
            && (is_extrapolation_rule_auto_constructible_v<
                        MaxOf<ExtrapRules>,
                        CoeffGridOf<Basis>,
                        DataType,
                        Basis> && ...))
        : m_evaluator(LagrangeEvaluator<
                      ExecSpace,
                      memory_space,
                      DataType,
                      Basis,
                      Grid,
                      MinOf<ExtrapRules>,
                      MaxOf<ExtrapRules>>(
                get_extrapolation<MinOf<ExtrapRules>, CoeffGridOf<Basis>, DataType, Basis>(
                        Extremity::FRONT),
                get_extrapolation<MaxOf<ExtrapRules>, CoeffGridOf<Basis>, DataType, Basis>(
                        Extremity::BACK))...)
    {
    }

    /**
     * @brief Construct an NDLagrangeInterpolator, specifying the extrapolation rules
     * explicitly.
     *
     * Use this overload when a chosen extrapolation rule cannot be built automatically,
     * e.g. a custom extrapolation rule that is not default-constructible and is not
     * ExtrapolationRule::Constant.
     *
     * @param idx_range The index range on which the interpolator will act. This is
     *                  unused but is included to match the SplineInterpolator interface.
     * @param extrapolation_rules One std::pair<MinExtrapolationRuleI, MaxExtrapolationRuleI>
     *                  per dimension, in the same order as IdxRangeBasis/IdxRangeInterpGrid.
     */
    explicit NDLagrangeInterpolator(
            IdxRange<Grid...> idx_range,
            std::pair<MinOf<ExtrapRules>, MaxOf<ExtrapRules>> const&... extrapolation_rules)
        : m_evaluator(LagrangeEvaluator<
                      ExecSpace,
                      memory_space,
                      DataType,
                      Basis,
                      Grid,
                      MinOf<ExtrapRules>,
                      MaxOf<ExtrapRules>>(extrapolation_rules.first, extrapolation_rules.second)...)
    {
    }

    /**
     * @brief Return a const reference to the owned ND identity builder.
     * @return The BuilderType instance.
     */
    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    /**
     * @brief Return a const reference to the owned ND Lagrange evaluator.
     * @return The EvaluatorType instance.
     */
    EvaluatorType const& get_evaluator() const
    {
        return m_evaluator;
    }
};

/**
 * @brief A helper alias to define an instance of detail::NDLagrangeInterpolator.
 *
 * The helper allows ExtrapRulesSeq to be more general. It is a
 * ddc::detail::TypeSeq<ExtrapRules1, ..., ExtrapRulesN> pairing, for each dimension, the
 * extrapolation rules applied below/above the boundary. Each ExtrapRulesI is itself a
 * ddc::detail::TypeSeq<MinExtrapolationRuleI, MaxExtrapolationRuleI> where each rule may
 * be one of the tags in the ExtrapolationRule namespace (e.g. ExtrapolationRule::Periodic)
 * or a custom, already-concrete extrapolation rule class.
 */
template <
        class ExecSpace,
        class DataType,
        class IdxRangeBasis,
        class IdxRangeInterpGrid,
        class... ExtrapRulesSeq>
struct LagrangeInterpolatorResolver;

template <class ExecSpace, class DataType, class Basis, class Grid, class ExtrapRules>
struct LagrangeInterpolatorResolver<
        ExecSpace,
        DataType,
        IdxRange<Basis>,
        IdxRange<Grid>,
        ExtrapRules>
{
    using type = detail::LagrangeInterpolator<
            ExecSpace,
            IdxRange<Basis>,
            IdxRange<Grid>,
            extrapolation_rule_t<
                    ExtrapRules,
                    typename IdentityInterpolationBuilder<
                            ExecSpace,
                            typename ExecSpace::memory_space,
                            DataType,
                            Grid,
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
        class GridHead,
        class... Grid,
        class ExtrapRulesHead,
        class... ExtrapRules>
struct LagrangeInterpolatorResolver<
        ExecSpace,
        DataType,
        IdxRange<BasisHead, Basis...>,
        IdxRange<GridHead, Grid...>,
        ExtrapRulesHead,
        ExtrapRules...>
{
    using type = detail::NDLagrangeInterpolator<
            ExecSpace,
            IdxRange<BasisHead, Basis...>,
            IdxRange<GridHead, Grid...>,
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

/**
 * @brief A helper alias to define an instance of detail::LagrangeInterpolator.
 *
 * The helper allows ExtrapRules to be more general. It is a
 * ddc::detail::TypeSeq<MinExtrapolationRule, MaxExtrapolationRule> pairing the
 * extrapolation rules applied below/above the boundary. Each may
 * be one of the tags in the ExtrapolationRule namespace (e.g.
 * ExtrapolationRule::Periodic) or a custom, already-concrete
 * extrapolation rule class.
 */
template <
        class ExecSpace,
        class DataType,
        class IdxRangeBasis,
        class IdxRangeInterpGrid,
        class... ExtrapRules>
using LagrangeInterpolator = typename detail::NDLagrangeInterpolatorResolver<
        ExecSpace,
        DataType,
        IdxRangeBasis,
        IdxRangeInterpGrid,
        ExtrapRules...>::type;
