

# File spline\_interpolation.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**spline\_interpolation.hpp**](spline__interpolation_8hpp.md)

[Go to the documentation of this file](spline__interpolation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <utility>

#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "extrapolation_rule_choice.hpp"

template <ddc::SplineBuilderClosure MinClosure, ddc::SplineBuilderClosure MaxClosure>
struct SplineBoundaryClosures
{
    using min = std::integral_constant<ddc::SplineBuilderClosure, MinClosure>;
    using max = std::integral_constant<ddc::SplineBuilderClosure, MaxClosure>;
};

namespace SplineBoundaryClosure {

using Periodic = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC>;

using Greville_Greville = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::GREVILLE,
        ddc::SplineBuilderClosure::GREVILLE>;

using Hermite_Hermite = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::HERMITE,
        ddc::SplineBuilderClosure::HERMITE>;

using HomogeneousHermite_HomogeneousHermite = SplineBoundaryClosures<
        ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE,
        ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE>;

} // namespace SplineBoundaryClosure

namespace detail {

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
    using BuilderType = ddc::SplineBuilder<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis,
            InterpGrid,
            MinBound,
            MaxBound,
            Solver>;

    using EvaluatorType = ddc::SplineEvaluator<
            ExecSpace,
            typename ExecSpace::memory_space,
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
        : m_min_extrapolation(get_extrapolation<
                              MinExtrapolationRule,
                              double,
                              typename Basis::continuous_dimension_type,
                              IdxRange<InterpGrid>,
                              IdxRange<Basis>>(Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<
                              MaxExtrapolationRule,
                              double,
                              typename Basis::continuous_dimension_type,
                              IdxRange<InterpGrid>,
                              IdxRange<Basis>>(Extremity::BACK))
        , m_builder(label, idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }
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
        : m_min_extrapolation(get_extrapolation<
                              MinExtrapolationRule,
                              double,
                              typename Basis::continuous_dimension_type,
                              IdxRange<InterpGrid>,
                              IdxRange<Basis>>(Extremity::FRONT))
        , m_max_extrapolation(get_extrapolation<
                              MaxExtrapolationRule,
                              double,
                              typename Basis::continuous_dimension_type,
                              IdxRange<InterpGrid>,
                              IdxRange<Basis>>(Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(m_min_extrapolation, m_max_extrapolation)
    {
    }

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
        class Basis1,
        class Basis2,
        class InterpGrid1,
        class InterpGrid2,
        class ExtrapRules,
        class BoundaryClosures,
        ddc::SplineSolver Solver>
class SplineInterpolator2D;

template <
        class ExecSpace,
        class Basis1,
        class Basis2,
        class InterpGrid1,
        class InterpGrid2,
        class ExtrapRules1,
        class ExtrapRules2,
        class BoundaryClosures1,
        class BoundaryClosures2,
        ddc::SplineSolver Solver>
class SplineInterpolator2D<
        ExecSpace,
        Basis1,
        Basis2,
        InterpGrid1,
        InterpGrid2,
        ddc::detail::TypeSeq<ExtrapRules1, ExtrapRules2>,
        ddc::detail::TypeSeq<BoundaryClosures1, BoundaryClosures2>,
        Solver>
{
private:
    using continuous_dimension_type1 = typename InterpGrid1::continuous_dimension_type;
    using continuous_dimension_type2 = typename InterpGrid2::continuous_dimension_type;

    static constexpr bool is_periodic1 = continuous_dimension_type1::PERIODIC;
    static constexpr bool is_periodic2 = continuous_dimension_type2::PERIODIC;

    using MinExtrapolationRule1 = ddc::type_seq_element_t<0, ExtrapRules1>;
    using MaxExtrapolationRule1 = ddc::type_seq_element_t<1, ExtrapRules1>;
    using MinExtrapolationRule2 = ddc::type_seq_element_t<0, ExtrapRules2>;
    using MaxExtrapolationRule2 = ddc::type_seq_element_t<1, ExtrapRules2>;

    static constexpr ddc::SplineBuilderClosure MinBound1 = BoundaryClosures1::min::value;
    static constexpr ddc::SplineBuilderClosure MaxBound1 = BoundaryClosures1::max::value;
    static constexpr ddc::SplineBuilderClosure MinBound2 = BoundaryClosures2::min::value;
    static constexpr ddc::SplineBuilderClosure MaxBound2 = BoundaryClosures2::max::value;

    static_assert(is_periodic1 == (MinBound1 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic1 == (MaxBound1 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic1
            == std::is_same_v<
                    MinExtrapolationRule1,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type1>>);
    static_assert(
            is_periodic1
            == std::is_same_v<
                    MaxExtrapolationRule1,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type1>>);

    static_assert(is_periodic2 == (MinBound2 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic2 == (MaxBound2 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic2
            == std::is_same_v<
                    MinExtrapolationRule2,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type2>>);
    static_assert(
            is_periodic2
            == std::is_same_v<
                    MaxExtrapolationRule2,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type2>>);

    static_assert(is_spline_basis_v<Basis1>);
    static_assert(is_spline_basis_v<Basis2>);

public:
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
    explicit SplineInterpolator2D(IdxRange<InterpGrid1, InterpGrid2> idx_range) requires(
            (is_extrapolation_rule_auto_constructible_v<MinExtrapolationRule1, InterpGrid1, double, Basis1>)&&(is_extrapolation_rule_auto_constructible_v<MaxExtrapolationRule1, InterpGrid1, double, Basis1>)&&(
                    is_extrapolation_rule_auto_constructible_v<
                            MinExtrapolationRule2,
                            InterpGrid2,
                            double,
                            Basis2>)&&(is_extrapolation_rule_auto_constructible_v<MaxExtrapolationRule2, InterpGrid2, double, Basis2>))
        : m_min_extrapolation1(get_extrapolation<
                               MinExtrapolationRule1,
                               double,
                               typename Basis1::continuous_dimension_type,
                               IdxRange<InterpGrid1, InterpGrid2>,
                               IdxRange<Basis1, Basis2>>(Extremity::FRONT))
        , m_max_extrapolation1(get_extrapolation<
                               MaxExtrapolationRule1,
                               double,
                               typename Basis1::continuous_dimension_type,
                               IdxRange<InterpGrid1, InterpGrid2>,
                               IdxRange<Basis1, Basis2>>(Extremity::BACK))
        , m_min_extrapolation2(get_extrapolation<
                               MinExtrapolationRule2,
                               double,
                               typename Basis2::continuous_dimension_type,
                               IdxRange<InterpGrid1, InterpGrid2>,
                               IdxRange<Basis1, Basis2>>(Extremity::FRONT))
        , m_max_extrapolation2(get_extrapolation<
                               MaxExtrapolationRule2,
                               double,
                               typename Basis2::continuous_dimension_type,
                               IdxRange<InterpGrid1, InterpGrid2>,
                               IdxRange<Basis1, Basis2>>(Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(
                  m_min_extrapolation1,
                  m_max_extrapolation1,
                  m_min_extrapolation2,
                  m_max_extrapolation2)
    {
    }

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
        class Basis1,
        class Basis2,
        class Basis3,
        class InterpGrid1,
        class InterpGrid2,
        class InterpGrid3,
        class ExtrapRules,
        class BoundaryClosures,
        ddc::SplineSolver Solver>
class SplineInterpolator3D;

template <
        class ExecSpace,
        class Basis1,
        class Basis2,
        class Basis3,
        class InterpGrid1,
        class InterpGrid2,
        class InterpGrid3,
        class ExtrapRules1,
        class ExtrapRules2,
        class ExtrapRules3,
        class BoundaryClosures1,
        class BoundaryClosures2,
        class BoundaryClosures3,
        ddc::SplineSolver Solver>
class SplineInterpolator3D<
        ExecSpace,
        Basis1,
        Basis2,
        Basis3,
        InterpGrid1,
        InterpGrid2,
        InterpGrid3,
        ddc::detail::TypeSeq<ExtrapRules1, ExtrapRules2, ExtrapRules3>,
        ddc::detail::TypeSeq<BoundaryClosures1, BoundaryClosures2, BoundaryClosures3>,
        Solver>
{
private:
    using continuous_dimension_type1 = typename InterpGrid1::continuous_dimension_type;
    using continuous_dimension_type2 = typename InterpGrid2::continuous_dimension_type;
    using continuous_dimension_type3 = typename InterpGrid3::continuous_dimension_type;

    static constexpr bool is_periodic1 = continuous_dimension_type1::PERIODIC;
    static constexpr bool is_periodic2 = continuous_dimension_type2::PERIODIC;
    static constexpr bool is_periodic3 = continuous_dimension_type3::PERIODIC;

    using MinExtrapolationRule1 = ddc::type_seq_element_t<0, ExtrapRules1>;
    using MaxExtrapolationRule1 = ddc::type_seq_element_t<1, ExtrapRules1>;
    using MinExtrapolationRule2 = ddc::type_seq_element_t<0, ExtrapRules2>;
    using MaxExtrapolationRule2 = ddc::type_seq_element_t<1, ExtrapRules2>;
    using MinExtrapolationRule3 = ddc::type_seq_element_t<0, ExtrapRules3>;
    using MaxExtrapolationRule3 = ddc::type_seq_element_t<1, ExtrapRules3>;

    static constexpr ddc::SplineBuilderClosure MinBound1 = BoundaryClosures1::min::value;
    static constexpr ddc::SplineBuilderClosure MaxBound1 = BoundaryClosures1::max::value;
    static constexpr ddc::SplineBuilderClosure MinBound2 = BoundaryClosures2::min::value;
    static constexpr ddc::SplineBuilderClosure MaxBound2 = BoundaryClosures2::max::value;
    static constexpr ddc::SplineBuilderClosure MinBound3 = BoundaryClosures3::min::value;
    static constexpr ddc::SplineBuilderClosure MaxBound3 = BoundaryClosures3::max::value;

    static_assert(is_periodic1 == (MinBound1 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic1 == (MaxBound1 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic1
            == std::is_same_v<
                    MinExtrapolationRule1,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type1>>);
    static_assert(
            is_periodic1
            == std::is_same_v<
                    MaxExtrapolationRule1,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type1>>);

    static_assert(is_periodic2 == (MinBound2 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic2 == (MaxBound2 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic2
            == std::is_same_v<
                    MinExtrapolationRule2,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type2>>);
    static_assert(
            is_periodic2
            == std::is_same_v<
                    MaxExtrapolationRule2,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type2>>);

    static_assert(is_periodic3 == (MinBound3 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(is_periodic3 == (MaxBound3 == ddc::SplineBuilderClosure::PERIODIC));
    static_assert(
            is_periodic3
            == std::is_same_v<
                    MinExtrapolationRule3,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type3>>);
    static_assert(
            is_periodic3
            == std::is_same_v<
                    MaxExtrapolationRule3,
                    ddc::PeriodicExtrapolationRule<continuous_dimension_type3>>);

    static_assert(is_spline_basis_v<Basis1>);
    static_assert(is_spline_basis_v<Basis2>);
    static_assert(is_spline_basis_v<Basis3>);

public:
    using BuilderType = ddc::SplineBuilder3D<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis1,
            Basis2,
            Basis3,
            InterpGrid1,
            InterpGrid2,
            InterpGrid3,
            MinBound1,
            MaxBound1,
            MinBound2,
            MaxBound2,
            MinBound3,
            MaxBound3,
            Solver>;

    using EvaluatorType = ddc::SplineEvaluator3D<
            ExecSpace,
            typename ExecSpace::memory_space,
            Basis1,
            Basis2,
            Basis3,
            InterpGrid1,
            InterpGrid2,
            InterpGrid3,
            MinExtrapolationRule1,
            MaxExtrapolationRule1,
            MinExtrapolationRule2,
            MaxExtrapolationRule2,
            MinExtrapolationRule3,
            MaxExtrapolationRule3>;

    static constexpr std::size_t rank()
    {
        return 3;
    }

private:
    MinExtrapolationRule1 m_min_extrapolation1;
    MaxExtrapolationRule1 m_max_extrapolation1;
    MinExtrapolationRule2 m_min_extrapolation2;
    MaxExtrapolationRule2 m_max_extrapolation2;
    MinExtrapolationRule3 m_min_extrapolation3;
    MaxExtrapolationRule3 m_max_extrapolation3;
    BuilderType m_builder;
    EvaluatorType m_evaluator;

public:
    explicit SplineInterpolator3D(IdxRange<
                                  InterpGrid1,
                                  InterpGrid2,
                                  InterpGrid3> idx_range) requires((is_extrapolation_rule_auto_constructible_v<MinExtrapolationRule1, InterpGrid1, double, Basis1>)&&(is_extrapolation_rule_auto_constructible_v<MaxExtrapolationRule1, InterpGrid1, double, Basis1>)&&(is_extrapolation_rule_auto_constructible_v<MinExtrapolationRule2, InterpGrid2, double, Basis2>)&&(is_extrapolation_rule_auto_constructible_v<MaxExtrapolationRule2, InterpGrid2, double, Basis2>)&&(is_extrapolation_rule_auto_constructible_v<MinExtrapolationRule3, InterpGrid3, double, Basis3>)&&(is_extrapolation_rule_auto_constructible_v<MaxExtrapolationRule3, InterpGrid3, double, Basis3>))
        : m_min_extrapolation1(
                get_extrapolation<MinExtrapolationRule1, InterpGrid1, double, Basis1>(
                        Extremity::FRONT))
        , m_max_extrapolation1(
                  get_extrapolation<MaxExtrapolationRule1, InterpGrid1, double, Basis1>(
                          Extremity::BACK))
        , m_min_extrapolation2(
                  get_extrapolation<MinExtrapolationRule2, InterpGrid2, double, Basis2>(
                          Extremity::FRONT))
        , m_max_extrapolation2(
                  get_extrapolation<MaxExtrapolationRule2, InterpGrid2, double, Basis2>(
                          Extremity::BACK))
        , m_min_extrapolation3(
                  get_extrapolation<MinExtrapolationRule3, InterpGrid3, double, Basis3>(
                          Extremity::FRONT))
        , m_max_extrapolation3(
                  get_extrapolation<MaxExtrapolationRule3, InterpGrid3, double, Basis3>(
                          Extremity::BACK))
        , m_builder(idx_range)
        , m_evaluator(
                  m_min_extrapolation1,
                  m_max_extrapolation1,
                  m_min_extrapolation2,
                  m_max_extrapolation2,
                  m_min_extrapolation3,
                  m_max_extrapolation3)
    {
    }

    explicit SplineInterpolator3D(
            IdxRange<InterpGrid1, InterpGrid2> idx_range,
            MinExtrapolationRule1 min_extrapolation_rule1,
            MaxExtrapolationRule1 max_extrapolation_rule1,
            MinExtrapolationRule2 min_extrapolation_rule2,
            MaxExtrapolationRule2 max_extrapolation_rule2,
            MinExtrapolationRule3 min_extrapolation_rule3,
            MaxExtrapolationRule3 max_extrapolation_rule3)
        : m_min_extrapolation1(std::move(min_extrapolation_rule1))
        , m_max_extrapolation1(std::move(max_extrapolation_rule1))
        , m_min_extrapolation2(std::move(min_extrapolation_rule2))
        , m_max_extrapolation2(std::move(max_extrapolation_rule2))
        , m_min_extrapolation3(std::move(min_extrapolation_rule3))
        , m_max_extrapolation3(std::move(max_extrapolation_rule3))
        , m_builder(idx_range)
        , m_evaluator(
                  m_min_extrapolation1,
                  m_max_extrapolation1,
                  m_min_extrapolation2,
                  m_max_extrapolation2,
                  m_min_extrapolation3,
                  m_max_extrapolation3)
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
        ddc::SplineSolver Solver,
        class... TailTags>
struct SplineInterpolatorResolver;

template <
        class ExecSpace,
        class Basis,
        class InterpGrid,
        ddc::SplineSolver Solver,
        class ExtrapRules,
        class BoundaryClosures>
struct SplineInterpolatorResolver<
        ExecSpace,
        IdxRange<Basis>,
        IdxRange<InterpGrid>,
        Solver,
        ExtrapRules,
        BoundaryClosures>
{
    using type = SplineInterpolator<
            ExecSpace,
            Basis,
            InterpGrid,
            extrapolation_rule_t<ExtrapRules, double, Basis, IdxRange<InterpGrid>>,
            BoundaryClosures,
            Solver>;
};

template <
        class ExecSpace,
        class Basis1,
        class Basis2,
        class InterpGrid1,
        class InterpGrid2,
        ddc::SplineSolver Solver,
        class ExtrapRules1,
        class ExtrapRules2,
        class BoundaryClosures1,
        class BoundaryClosures2>
struct SplineInterpolatorResolver<
        ExecSpace,
        IdxRange<Basis1, Basis2>,
        IdxRange<InterpGrid1, InterpGrid2>,
        Solver,
        ExtrapRules1,
        ExtrapRules2,
        BoundaryClosures1,
        BoundaryClosures2>
{
    using type = SplineInterpolator2D<
            ExecSpace,
            Basis1,
            Basis2,
            InterpGrid1,
            InterpGrid2,
            ddc::detail::TypeSeq<
                    extrapolation_rule_t<
                            ExtrapRules1,
                            double,
                            Basis1,
                            IdxRange<InterpGrid1, InterpGrid2>>,
                    extrapolation_rule_t<
                            ExtrapRules2,
                            double,
                            Basis2,
                            IdxRange<InterpGrid1, InterpGrid2>>>,
            ddc::detail::TypeSeq<BoundaryClosures1, BoundaryClosures2>,
            Solver>;
};

template <
        class ExecSpace,
        class Basis1,
        class Basis2,
        class Basis3,
        class InterpGrid1,
        class InterpGrid2,
        class InterpGrid3,
        ddc::SplineSolver Solver,
        class ExtrapRules1,
        class ExtrapRules2,
        class ExtrapRules3,
        class BoundaryClosures1,
        class BoundaryClosures2,
        class BoundaryClosures3>
struct SplineInterpolatorResolver<
        ExecSpace,
        IdxRange<Basis1, Basis2, Basis3>,
        IdxRange<InterpGrid1, InterpGrid2, InterpGrid3>,
        Solver,
        ExtrapRules1,
        ExtrapRules2,
        ExtrapRules3,
        BoundaryClosures1,
        BoundaryClosures2,
        BoundaryClosures3>
{
    using type = SplineInterpolator3D<
            ExecSpace,
            Basis1,
            Basis2,
            Basis3,
            InterpGrid1,
            InterpGrid2,
            InterpGrid3,
            ddc::detail::TypeSeq<
                    extrapolation_rule_t<
                            ExtrapRules1,
                            double,
                            Basis1,
                            IdxRange<InterpGrid1, InterpGrid2, InterpGrid3>>,
                    extrapolation_rule_t<
                            ExtrapRules2,
                            double,
                            Basis2,
                            IdxRange<InterpGrid1, InterpGrid2, InterpGrid3>>,
                    extrapolation_rule_t<
                            ExtrapRules3,
                            double,
                            Basis3,
                            IdxRange<InterpGrid1, InterpGrid2, InterpGrid3>>>,
            ddc::detail::TypeSeq<BoundaryClosures1, BoundaryClosures2, BoundaryClosures3>,
            Solver>;
};

} // namespace detail

template <class ExecSpace, class IdxRangeBasis, class IdxRangeInterpGrid, class... TailTags>
using SplineInterpolator = typename detail::SplineInterpolatorResolver<
        ExecSpace,
        IdxRangeBasis,
        IdxRangeInterpGrid,
        ddc::SplineSolver::LAPACK,
        TailTags...>::type;

template <
        class ExecSpace,
        class IdxRangeBasis,
        class IdxRangeInterpGrid,
        ddc::SplineSolver Solver,
        class... TailTags>
using SplineInterpolatorWSolver = typename detail::SplineInterpolatorResolver<
        ExecSpace,
        IdxRangeBasis,
        IdxRangeInterpGrid,
        Solver,
        TailTags...>::type;
```


