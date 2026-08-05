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
 * @tparam MinExtrapRule1 The extrapolation rule applied below the lower boundary of the
 *                        first dimension. This may be one of the tags in the
 *                        ExtrapolationRule namespace (e.g. ExtrapolationRule::Periodic)
 *                        or a custom, already-concrete extrapolation rule class.
 * @tparam MaxExtrapRule1 The extrapolation rule applied above the upper boundary of the
 *                        first dimension. See MinExtrapRule1 for the accepted forms.
 * @tparam MinExtrapRule2 The extrapolation rule applied below the lower boundary of the
 *                        second dimension. See MinExtrapRule1 for the accepted forms.
 * @tparam MaxExtrapRule2 The extrapolation rule applied above the upper boundary of the
 *                        second dimension. See MinExtrapRule1 for the accepted forms.
 * @tparam MinBound1      The ddc::SplineBuilderClosure at the lower boundary of the first dimension.
 * @tparam MaxBound1      The ddc::SplineBuilderClosure at the upper boundary of the first dimension.
 * @tparam MinBound2      The ddc::SplineBuilderClosure at the lower boundary of the second dimension.
 * @tparam MaxBound2      The ddc::SplineBuilderClosure at the upper boundary of the second dimension.
 * @tparam Solver         The spline solver backend (default: LAPACK).
 */
template <
        class ExecSpace,
        class Basis1,
        class Basis2,
        class InterpGrid1,
        class InterpGrid2,
        class MinExtrapRule1,
        class MaxExtrapRule1,
        class MinExtrapRule2,
        class MaxExtrapRule2,
        ddc::SplineBuilderClosure MinBound1,
        ddc::SplineBuilderClosure MaxBound1,
        ddc::SplineBuilderClosure MinBound2,
        ddc::SplineBuilderClosure MaxBound2,
        ddc::SplineSolver Solver = ddc::SplineSolver::LAPACK>
class SplineInterpolator2D
{
private:
    using continuous_dimension_type1 = typename InterpGrid1::continuous_dimension_type;
    using continuous_dimension_type2 = typename InterpGrid2::continuous_dimension_type;

    static constexpr bool is_periodic1 = continuous_dimension_type1::PERIODIC;
    static constexpr bool is_periodic2 = continuous_dimension_type2::PERIODIC;

    static_assert(is_periodic1 == (MinBound1 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic1 == (MaxBound1 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic1
            == std::is_same_v<
                    extrapolation_rule_t<MinExtrapRule1, Basis1, double>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type1>>);
    static_assert(
            is_periodic1
            == std::is_same_v<
                    extrapolation_rule_t<MaxExtrapRule1, Basis1, double>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type1>>);

    static_assert(is_periodic2 == (MinBound2 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic2 == (MaxBound2 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic2
            == std::is_same_v<
                    extrapolation_rule_t<MinExtrapRule2, Basis2, double>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type2>>);
    static_assert(
            is_periodic2
            == std::is_same_v<
                    extrapolation_rule_t<MaxExtrapRule2, Basis2, double>,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type2>>);

    static_assert(is_spline_basis_v<Basis1>);
    static_assert(is_spline_basis_v<Basis2>);

    using MinExtrapolationRule1 = extrapolation_rule_t<MinExtrapRule1, Basis1, double>;
    using MaxExtrapolationRule1 = extrapolation_rule_t<MaxExtrapRule1, Basis1, double>;
    using MinExtrapolationRule2 = extrapolation_rule_t<MinExtrapRule2, Basis2, double>;
    using MaxExtrapolationRule2 = extrapolation_rule_t<MaxExtrapRule2, Basis2, double>;

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
            MinExtrapolationRule1,
            MaxExtrapolationRule1,
            MinExtrapolationRule2,
            MaxExtrapolationRule2>;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return 2;
    }

private:
    MinExtrapolationRule1 m_min_extrapolation1;
    MaxExtrapolationRule1 m_max_extrapolation1;
    MinExtrapolationRule2 m_min_extrapolation2;
    MaxExtrapolationRule2 m_max_extrapolation2;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    /**
     * @brief Construct a SplineInterpolator2D on the given 2D interpolation index range.
     *
     * The extrapolation rules are initialised from the discrete spaces of @c Basis1 and
     * @c Basis2, so the corresponding ddc discrete spaces must be initialised before
     * construction. This overload is only available when all four extrapolation rules
     * can be built automatically (they are default-constructible, or the tag
     * ExtrapolationRule::Constant is used) - otherwise use the overload that takes the
     * extrapolation rules explicitly.
     *
     * @param idx_range The 2D interpolation index range passed to the builder.
     */
    explicit SplineInterpolator2D(IdxRange<InterpGrid1, InterpGrid2> idx_range) requires(
            is_extrapolation_rule_auto_constructible_v<MinExtrapRule1, Basis1, double>&&
                    is_extrapolation_rule_auto_constructible_v<MaxExtrapRule1, Basis1, double>&&
                            is_extrapolation_rule_auto_constructible_v<
                                    MinExtrapRule2,
                                    Basis2,
                                    double>&&
                                    is_extrapolation_rule_auto_constructible_v<
                                            MaxExtrapRule2,
                                            Basis2,
                                            double>)
        : m_min_extrapolation1(get_extrapolation<MinExtrapRule1, Basis1, double>(Extremity::FRONT))
        , m_max_extrapolation1(get_extrapolation<MaxExtrapRule1, Basis1, double>(Extremity::BACK))
        , m_min_extrapolation2(get_extrapolation<MinExtrapRule2, Basis2, double>(Extremity::FRONT))
        , m_max_extrapolation2(get_extrapolation<MaxExtrapRule2, Basis2, double>(Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(
                  m_min_extrapolation1,
                  m_max_extrapolation1,
                  m_min_extrapolation2,
                  m_max_extrapolation2)
    {
    }

    /**
     * @brief Construct a SplineInterpolator2D on the given 2D interpolation index range,
     * specifying the extrapolation rules explicitly.
     *
     * Use this overload when a chosen extrapolation rule cannot be built automatically,
     * e.g. a custom extrapolation rule that is not default-constructible and is not
     * ExtrapolationRule::Constant.
     *
     * @param idx_range The 2D interpolation index range passed to the builder.
     * @param min_extrapolation_rule1 The extrapolation rule to use below the lower
     *                                boundary of the first dimension.
     * @param max_extrapolation_rule1 The extrapolation rule to use above the upper
     *                                boundary of the first dimension.
     * @param min_extrapolation_rule2 The extrapolation rule to use below the lower
     *                                boundary of the second dimension.
     * @param max_extrapolation_rule2 The extrapolation rule to use above the upper
     *                                boundary of the second dimension.
     */
    explicit SplineInterpolator2D(
            IdxRange<InterpGrid1, InterpGrid2> idx_range,
            MinExtrapolationRule1 min_extrapolation_rule1,
            MaxExtrapolationRule1 max_extrapolation_rule1,
            MinExtrapolationRule2 min_extrapolation_rule2,
            MaxExtrapolationRule2 max_extrapolation_rule2)
        : m_min_extrapolation1(std::move(min_extrapolation_rule1))
        , m_max_extrapolation1(std::move(max_extrapolation_rule1))
        , m_min_extrapolation2(std::move(min_extrapolation_rule2))
        , m_max_extrapolation2(std::move(max_extrapolation_rule2))
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
