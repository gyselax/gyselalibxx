// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "extrapolation_rule_choice.hpp"

/**
 * @brief An owning interpolation object that bundles a spline builder and evaluator.
 *
 * SplineInterpolator constructs and owns a matching ddc::SplineBuilder and
 * ddc::SplineEvaluator for a given dimension. It satisfies the
 * concepts::Interpolation concept and is the recommended way to create a
 * spline interpolation for use with advection operators and similar algorithms.
 *
 * The boundary condition (MinBound / MaxBound) and extrapolation rule
 * (MinExtrapRule / MaxExtrapRule) must be consistent: both must be PERIODIC for
 * periodic dimensions and both must be non-PERIODIC for non-periodic dimensions.
 *
 * @tparam ExecSpace     The Kokkos execution space used for computations.
 * @tparam Basis         The B-spline basis type (uniform or non-uniform).
 * @tparam InterpGrid    The discrete grid on which function values are provided.
 * @tparam MinExtrapRule The ExtrapolationRule applied below the lower boundary.
 * @tparam MaxExtrapRule The ExtrapolationRule applied above the upper boundary.
 * @tparam MinBound      The ddc::BoundCond at the lower boundary of the spline builder.
 * @tparam MaxBound      The ddc::BoundCond at the upper boundary of the spline builder.
 * @tparam Solver        The spline solver backend (default: LAPACK).
 */
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
    /// @brief The ddc::SplineBuilder type built from the template parameters.
    using BuilderType = ddc::SplineBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis,
            InterpGrid,
            MinBound,
            MaxBound,
            Solver>;

    /// @brief The ddc::SplineEvaluator type built from the template parameters.
    using EvaluatorType = ddc::SplineEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis,
            InterpGrid,
            extrapolation_rule_t<MinExtrapRule, Basis>,
            extrapolation_rule_t<MaxExtrapRule, Basis>>;

    /// @brief The number of interpolation dimensions.
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
    /**
     * @brief Construct a SplineInterpolator on the given interpolation index range.
     *
     * The extrapolation rules are initialised from the discrete space of @c Basis,
     * so the corresponding ddc discrete space must be initialised before construction.
     *
     * @param idx_range The 1D interpolation index range passed to the builder.
     */
    explicit SplineInterpolator(IdxRange<InterpGrid> idx_range)
        : m_min_extrapolation(get_extrapolation<MinExtrapRule, Basis>(Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<MaxExtrapRule, Basis>(Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

    /**
     * @brief Return a const reference to the owned spline builder.
     * @return The BuilderType instance.
     */
    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    /**
     * @brief Return a const reference to the owned spline evaluator.
     * @return The EvaluatorType instance.
     */
    EvaluatorType const& get_evaluator() const
    {
        return m_evaluator;
    }
};
