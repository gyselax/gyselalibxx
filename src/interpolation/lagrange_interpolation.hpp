// SPDX-License-Identifier: MIT
#pragma once

#include <utility>

#include "extrapolation_rule_choice.hpp"
#include "identity_interpolation_builder.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "lagrange_evaluator.hpp"

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
        class Basis,
        class InterpGrid,
        class ExtrapRules,
        class DataType = double>
using LagrangeInterpolator = detail::LagrangeInterpolator<
        ExecSpace,
        Basis,
        InterpGrid,
        extrapolation_rule_t<
                ExtrapRules,
                typename IdentityInterpolationBuilder<
                        ExecSpace,
                        typename ExecSpace::memory_space,
                        DataType,
                        InterpGrid,
                        Basis>::basis_domain_type,
                DataType,
                Basis>,
        DataType>;
