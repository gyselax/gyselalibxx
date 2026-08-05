// SPDX-License-Identifier: MIT
#pragma once
#include <utility>

#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "extrapolation_rule_choice.hpp"

/**
 * @brief Groups the lower (min) and upper (max) ddc::SplineBuilderClosure of a
 * spline builder into a single non-type template argument.
 */
template <ddc::SplineBuilderClosure MinClosure, ddc::SplineBuilderClosure MaxClosure>
struct SplineBoundaryClosures
{
    /// @brief The ddc::SplineBuilderClosure at the lower boundary of the spline builder.
    using min = std::integral_constant<ddc::SplineBuilderClosure, MinClosure>;
    /// @brief The ddc::SplineBuilderClosure at the upper boundary of the spline builder.
    using max = std::integral_constant<ddc::SplineBuilderClosure, MaxClosure>;
};

/// @brief Predefined SplineBoundaryClosures for the common case where the same
/// closure applies at both boundaries.
namespace SplineBoundaryClosure {

/// @brief Convenience pairing of ddc::SplineBuilderClosure::PERIODIC for both boundaries.
using Periodic = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC>;

/// @brief Convenience pairing of ddc::SplineBuilderClosure::GREVILLE for both boundaries.
using Greville_Greville = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::GREVILLE,
        ddc::SplineBuilderClosure::GREVILLE>;

/// @brief Convenience pairing of ddc::SplineBuilderClosure::HERMITE for both boundaries.
using Hermite_Hermite = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::HERMITE>;

/// @brief Convenience pairing of ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE for both
/// boundaries.
using HomogeneousHermite_HomogeneousHermite = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE,
        ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE>;

} // namespace SplineBoundaryClosure

namespace detail {

/**
 * @brief An owning interpolation object that bundles a spline builder and evaluator.
 *
 * SplineInterpolator constructs and owns a matching ddc::SplineBuilder and
 * ddc::SplineEvaluator for a given dimension. It satisfies the
 * concepts::Interpolation concept and is the recommended way to create a
 * spline interpolation for use with advection operators and similar algorithms.
 *
 * The boundary condition (BoundaryClosures) and extrapolation rule
 * (the Min/Max pair in ExtrapRules) must be consistent: both must be PERIODIC for
 * periodic dimensions and both must be non-PERIODIC for non-periodic dimensions.
 *
 * @tparam ExecSpace     The Kokkos execution space used for computations.
 * @tparam Basis         The B-spline basis type (uniform or non-uniform).
 * @tparam InterpGrid    The discrete grid on which function values are provided.
 * @tparam ExtrapRules   A ddc::detail::TypeSeq<MinExtrapolationRule, MaxExtrapolationRule> pairing the
 *                       extrapolation rules applied below/above the boundary. Where
 *                       MinExtrapolationRule and MaxExtrapolationRule are extrapolation rule classes.
 * @tparam BoundaryClosures A SplineBoundaryClosures pairing the ddc::SplineBuilderClosure
 *                       applied at the lower/upper boundary of the spline builder.
 * @tparam Solver        The spline solver backend (default: LAPACK).
 */
template <
        class ExecSpace,
        class Basis,
        class InterpGrid,
        class ExtrapRules,
        class BoundaryClosures,
        ddc::SplineSolver Solver>
class SplineInterpolator
{
private:
    using continuous_dimension_type = typename InterpGrid::continuous_dimension_type;

    using MinExtrapolationRule = ddc::type_seq_element_t<0, ExtrapRules>;
    using MaxExtrapolationRule = ddc::type_seq_element_t<1, ExtrapRules>;

    static constexpr ddc::SplineBuilderClosure MinBound = BoundaryClosures::min::value;
    static constexpr ddc::SplineBuilderClosure MaxBound = BoundaryClosures::max::value;

    static constexpr bool is_periodic = continuous_dimension_type::PERIODIC;

    static_assert(is_periodic == (MinBound == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic == (MaxBound == ddc::SplineBuilderClosure::PERIODIC));
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
     * @brief Construct a SplineInterpolator on the given interpolation index range.
     *
     * The extrapolation rules are initialised from the discrete space of @c Basis,
     * so the corresponding ddc discrete space must be initialised before construction.
     * This overload is only available when both extrapolation rules can be built
     * automatically (they are default-constructible, or the tag ExtrapolationRule::Constant
     * is used) - otherwise use the overload that takes the extrapolation rules explicitly.
     *
     * @param label A label used to tag parallel regions and memory allocations for profiling.
     * @param idx_range The 1D interpolation index range passed to the builder.
     */
    explicit SplineInterpolator(std::string const& label, IdxRange<InterpGrid> idx_range) requires(
            is_extrapolation_rule_auto_constructible_v<
                    MinExtrapolationRule,
                    InterpGrid,
                    double,
                    Basis>&&
                    is_extrapolation_rule_auto_constructible_v<
                            MaxExtrapolationRule,
                            InterpGrid,
                            double,
                            Basis>)
        : m_min_extrapolation(get_extrapolation<MinExtrapolationRule, InterpGrid, double, Basis>(
                Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<MaxExtrapolationRule, InterpGrid, double, Basis>(
                  Extremity::BACK))
        , m_builder(label, idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }
    /**
     * @brief Construct a SplineInterpolator on the given interpolation index range.
     *
     * The extrapolation rules are initialised from the discrete space of @c Basis,
     * so the corresponding ddc discrete space must be initialised before construction.
     * This overload is only available when both extrapolation rules can be built
     * automatically (they are default-constructible, or the tag ExtrapolationRule::Constant
     * is used) - otherwise use the overload that takes the extrapolation rules explicitly.
     *
     * @param idx_range The 1D interpolation index range passed to the builder.
     */
    explicit SplineInterpolator(IdxRange<InterpGrid> idx_range) requires(
            is_extrapolation_rule_auto_constructible_v<
                    MinExtrapolationRule,
                    InterpGrid,
                    double,
                    Basis>&&
                    is_extrapolation_rule_auto_constructible_v<
                            MaxExtrapolationRule,
                            InterpGrid,
                            double,
                            Basis>)
        : m_min_extrapolation(get_extrapolation<MinExtrapolationRule, InterpGrid, double, Basis>(
                Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<MaxExtrapolationRule, InterpGrid, double, Basis>(
                  Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

    /**
     * @brief Construct a SplineInterpolator on the given interpolation index range,
     * specifying the extrapolation rules explicitly.
     *
     * Use this overload when the chosen extrapolation rule cannot be built
     * automatically, e.g. a custom extrapolation rule that is not default-constructible
     * and is not ExtrapolationRule::Constant.
     *
     * @param label A label used to tag parallel regions and memory allocations for profiling.
     * @param idx_range The 1D interpolation index range passed to the builder.
     * @param min_extrapolation_rule The extrapolation rule to use below the lower boundary.
     * @param max_extrapolation_rule The extrapolation rule to use above the upper boundary.
     */
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

    /**
     * @brief Construct a SplineInterpolator on the given interpolation index range,
     * specifying the extrapolation rules explicitly.
     *
     * Use this overload when the chosen extrapolation rule cannot be built
     * automatically, e.g. a custom extrapolation rule that is not default-constructible
     * and is not ExtrapolationRule::Constant.
     *
     * @param idx_range The 1D interpolation index range passed to the builder.
     * @param min_extrapolation_rule The extrapolation rule to use below the lower boundary.
     * @param max_extrapolation_rule The extrapolation rule to use above the upper boundary.
     */
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

} // namespace detail

/**
 * @brief A helper alias to define an instance of detail::SplineInterpolator.
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
        class BoundaryClosures,
        ddc::SplineSolver Solver = ddc::SplineSolver::LAPACK>
using SplineInterpolator = detail::SplineInterpolator<
        ExecSpace,
        Basis,
        InterpGrid,
        extrapolation_rule_t<ExtrapRules, InterpGrid, double, Basis>,
        BoundaryClosures,
        Solver>;
