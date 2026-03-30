// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "extrapolation_rule_choice.hpp"

/**
 * @brief An owning interpolation object that bundles a 2D spline builder and evaluator.
 *
 * SplineInterpolator2D constructs and owns a matching ddc::SplineBuilder2D and
 * ddc::SplineEvaluator2D for two given dimensions. It is the recommended way to
 * create a 2D spline interpolation for use with advection operators and similar
 * algorithms.
 *
 * The boundary conditions (MinBound1/MaxBound1 and MinBound2/MaxBound2) and
 * extrapolation rules (MinExtrapRule1/MaxExtrapRule1 and MinExtrapRule2/MaxExtrapRule2)
 * must be consistent within each dimension: both must be PERIODIC for periodic
 * dimensions and both must be non-PERIODIC for non-periodic dimensions.
 *
 * @tparam ExecSpace      The Kokkos execution space used for computations.
 * @tparam Basis1         The B-spline basis type for the first dimension (uniform or non-uniform).
 * @tparam Basis2         The B-spline basis type for the second dimension (uniform or non-uniform).
 * @tparam InterpGrid1    The discrete grid on which function values are provided along the first dimension.
 * @tparam InterpGrid2    The discrete grid on which function values are provided along the second dimension.
 * @tparam MinExtrapRule1 The ExtrapolationRule applied below the lower boundary of the first dimension.
 * @tparam MaxExtrapRule1 The ExtrapolationRule applied above the upper boundary of the first dimension.
 * @tparam MinExtrapRule2 The ExtrapolationRule applied below the lower boundary of the second dimension.
 * @tparam MaxExtrapRule2 The ExtrapolationRule applied above the upper boundary of the second dimension.
 * @tparam MinBound1      The ddc::BoundCond at the lower boundary of the first dimension.
 * @tparam MaxBound1      The ddc::BoundCond at the upper boundary of the first dimension.
 * @tparam MinBound2      The ddc::BoundCond at the lower boundary of the second dimension.
 * @tparam MaxBound2      The ddc::BoundCond at the upper boundary of the second dimension.
 * @tparam Solver         The spline solver backend (default: LAPACK).
 */
template <
        class ExecSpace,
        class Basis1,
        class Basis2,
        class InterpGrid1,
        class InterpGrid2,
        ExtrapolationRule MinExtrapRule1,
        ExtrapolationRule MaxExtrapRule1,
        ExtrapolationRule MinExtrapRule2,
        ExtrapolationRule MaxExtrapRule2,
        ddc::BoundCond MinBound1,
        ddc::BoundCond MaxBound1,
        ddc::BoundCond MinBound2,
        ddc::BoundCond MaxBound2,
        ddc::SplineSolver Solver = ddc::SplineSolver::LAPACK>
class SplineInterpolator2D
{
private:
    using continuous_dimension_type1 = typename InterpGrid1::continuous_dimension_type;
    using continuous_dimension_type2 = typename InterpGrid2::continuous_dimension_type;

    static constexpr bool is_periodic1 = continuous_dimension_type1::PERIODIC;
    static constexpr bool is_periodic2 = continuous_dimension_type2::PERIODIC;

    static_assert(is_periodic1 == (MinBound1 == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic1 == (MaxBound1 == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic1 == (MinExtrapRule1 == ExtrapolationRule::PERIODIC));
    static_assert(is_periodic1 == (MaxExtrapRule1 == ExtrapolationRule::PERIODIC));

    static_assert(is_periodic2 == (MinBound2 == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic2 == (MaxBound2 == ddc::BoundCond::PERIODIC));
    static_assert(is_periodic2 == (MinExtrapRule2 == ExtrapolationRule::PERIODIC));
    static_assert(is_periodic2 == (MaxExtrapRule2 == ExtrapolationRule::PERIODIC));

    static_assert(is_spline_basis_v<Basis1>);
    static_assert(is_spline_basis_v<Basis2>);

public:
    /// @brief The ddc::SplineBuilder2D type built from the template parameters.
    using BuilderType = ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis1,
            Basis2,
            InterpGrid1,
            InterpGrid2,
            MinBound1,
            MaxBound1,
            MinBound2,
            MaxBound2,
            Solver>;

    /// @brief The ddc::SplineEvaluator2D type built from the template parameters.
    using EvaluatorType = ddc::SplineEvaluator2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis1,
            Basis2,
            InterpGrid1,
            InterpGrid2,
            extrapolation_rule_t<MinExtrapRule1, Basis1>,
            extrapolation_rule_t<MaxExtrapRule1, Basis1>,
            extrapolation_rule_t<MinExtrapRule2, Basis2>,
            extrapolation_rule_t<MaxExtrapRule2, Basis2>>;

private:
    extrapolation_rule_t<MinExtrapRule1, Basis1> m_min_extrapolation1;
    extrapolation_rule_t<MaxExtrapRule1, Basis1> m_max_extrapolation1;
    extrapolation_rule_t<MinExtrapRule2, Basis2> m_min_extrapolation2;
    extrapolation_rule_t<MaxExtrapRule2, Basis2> m_max_extrapolation2;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    /**
     * @brief Construct a SplineInterpolator2D on the given 2D interpolation index range.
     *
     * The extrapolation rules are initialised from the discrete spaces of @c Basis1 and
     * @c Basis2, so the corresponding ddc discrete spaces must be initialised before
     * construction.
     *
     * @param idx_range The 2D interpolation index range passed to the builder.
     */
    explicit SplineInterpolator2D(IdxRange<InterpGrid1, InterpGrid2> idx_range)
        : m_min_extrapolation1(get_extrapolation<MinExtrapRule1, Basis1>(Extremity::FRONT))
        , m_max_extrapolation1(get_extrapolation<MaxExtrapRule1, Basis1>(Extremity::BACK))
        , m_min_extrapolation2(get_extrapolation<MinExtrapRule2, Basis2>(Extremity::FRONT))
        , m_max_extrapolation2(get_extrapolation<MaxExtrapRule2, Basis2>(Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(
                  m_min_extrapolation1,
                  m_max_extrapolation1,
                  m_min_extrapolation2,
                  m_max_extrapolation2)
    {
    }

    /**
     * @brief Return a const reference to the owned 2D spline builder.
     * @return The BuilderType instance.
     */
    BuilderType const& get_builder() const
    {
        return m_builder;
    }

    /**
     * @brief Return a const reference to the owned 2D spline evaluator.
     * @return The EvaluatorType instance.
     */
    EvaluatorType const& get_evaluator() const
    {
        return m_evaluator;
    }
};
