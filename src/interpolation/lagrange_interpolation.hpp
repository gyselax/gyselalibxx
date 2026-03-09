// SPDX-License-Identifier: MIT
#pragma once

#include "extrapolation_rule_choice.hpp"
#include "identity_interpolation_builder.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "lagrange_evaluator.hpp"

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
 * The boundary condition (MinBound / MaxBound) and extrapolation rule
 * (MinExtrapRule / MaxExtrapRule) must be consistent: both must be PERIODIC for
 * periodic dimensions and both must be non-PERIODIC for non-periodic dimensions.
 * Note: @c CONSTANT extrapolation is not supported for Lagrange interpolation.
 *
 * @tparam ExecSpace     The Kokkos execution space used for computations.
 * @tparam Basis         The Lagrange basis type (uniform or non-uniform).
 * @tparam InterpGrid    The discrete grid on which function values are provided.
 * @tparam MinExtrapRule The ExtrapolationRule applied below the lower boundary.
 * @tparam MaxExtrapRule The ExtrapolationRule applied above the upper boundary.
 * @tparam MinBound      The ddc::BoundCond at the lower boundary (default: GREVILLE).
 *                       This is included to have an interface interchangeable with SplineBuilder
 *                       but is unused.
 * @tparam MaxBound      The ddc::BoundCond at the upper boundary (default: GREVILLE).
 *                       This is included to have an interface interchangeable with SplineBuilder
 *                       but is unused.
 * @tparam DataType      The floating-point type of the function values (default: double).
 */
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
    /// @brief The IdentityInterpolationBuilder type built from the template parameters.
    using BuilderType = IdentityInterpolationBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            DataType,
            InterpGrid,
            Basis>;

    /// @brief The discrete grid type used for the Lagrange coefficients (the Lagrange basis grid).
    using CoeffGridType = typename BuilderType::basis_domain_type;

    /// @brief The LagrangeEvaluator type built from the template parameters.
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
    /**
     * @brief Construct a LagrangeInterpolator.
     *
     * The extrapolation rules are initialised from the discrete space of @c Basis,
     * so the corresponding ddc discrete space must be initialised before construction.
     * No index range is required because the identity builder needs none.
     */
    LagrangeInterpolator()
        : m_min_extrapolation(
                get_extrapolation<MinExtrapRule, CoeffGridType, Basis>(Extremity::FRONT))
        , m_max_extrapolation(
                  get_extrapolation<MaxExtrapRule, CoeffGridType, Basis>(Extremity::BACK))
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
