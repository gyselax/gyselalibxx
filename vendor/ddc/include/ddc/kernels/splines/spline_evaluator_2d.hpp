// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include "Kokkos_Macros.hpp"
#include "periodic_extrapolation_rule.hpp"
#include "spline_boundary_conditions.hpp"
#include "view.hpp"

namespace ddc {

/**
 * @brief A class to evaluate, differentiate or integrate a 2D spline function.
 *
 * A class which contains an operator () which can be used to evaluate, differentiate or integrate a 2D spline function.
 *
 * @tparam ExecSpace The Kokkos execution space on which the spline evaluation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data (spline coefficients and evaluation) is stored.
 * @tparam BSplines1 The discrete dimension representing the B-splines along the first dimension of interest.
 * @tparam BSplines2 The discrete dimension representing the B-splines along the second dimension of interest.
 * @tparam EvaluationMesh1 The first discrete dimension on which evaluation points are defined.
 * @tparam EvaluationMesh2 The second discrete dimension on which evaluation points are defined.
 * @tparam LeftExtrapolationRule1 The lower extrapolation rule type along first dimension of interest.
 * @tparam RightExtrapolationRule1 The upper extrapolation rule type along first dimension of interest.
 * @tparam LeftExtrapolationRule2 The lower extrapolation rule type along second dimension of interest.
 * @tparam RightExtrapolationRule2 The upper extrapolation rule type along second dimension of interest.
 * @tparam IDimX A variadic template of all the discrete dimensions forming the full space (EvaluationMesh1 + EvaluationMesh2 + batched dimensions).
 */
template <
        class ExecSpace,
        class MemorySpace,
        class BSplines1,
        class BSplines2,
        class EvaluationMesh1,
        class EvaluationMesh2,
        class LeftExtrapolationRule1,
        class RightExtrapolationRule1,
        class LeftExtrapolationRule2,
        class RightExtrapolationRule2,
        class... IDimX>
class SplineEvaluator2D
{
private:
    /**
     * @brief Tag to indicate that the value of the spline should be evaluated.
     */
    struct eval_type
    {
    };

    /**
     * @brief Tag to indicate that derivative of the spline should be evaluated.
     */
    struct eval_deriv_type
    {
    };

    using tag_type1 = typename BSplines1::tag_type;
    using tag_type2 = typename BSplines2::tag_type;

public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;

    /// @brief The type of the first discrete dimension of interest used by this class.
    using evaluation_mesh_type1 = EvaluationMesh1;

    /// @brief The type of the second discrete dimension of interest used by this class.
    using evaluation_mesh_type2 = EvaluationMesh2;

    /// @brief The discrete dimension representing the B-splines along first dimension.
    using bsplines_type1 = BSplines1;

    /// @brief The discrete dimension representing the B-splines along second dimension.
    using bsplines_type2 = BSplines2;

    /// @brief The type of the domain for the 1D evaluation mesh along first dimension used by this class.
    using evaluation_domain_type1 = ddc::DiscreteDomain<evaluation_mesh_type1>;

    /// @brief The type of the domain for the 1D evaluation mesh along second dimension used by this class.
    using evaluation_domain_type2 = ddc::DiscreteDomain<evaluation_mesh_type2>;

    /// @brief The type of the domain for the 2D evaluation mesh used by this class.
    using evaluation_domain_type
            = ddc::DiscreteDomain<evaluation_mesh_type1, evaluation_mesh_type2>;

    /// @brief The type of the whole domain representing evaluation points.
    using batched_evaluation_domain_type = ddc::DiscreteDomain<IDimX...>;

    /// @brief The type of the 1D spline domain corresponding to the first dimension of interest.
    using spline_domain_type1 = ddc::DiscreteDomain<bsplines_type1>;

    /// @brief The type of the 1D spline domain corresponding to the second dimension of interest.
    using spline_domain_type2 = ddc::DiscreteDomain<bsplines_type2>;

    /// @brief The type of the 2D spline domain corresponding to the dimensions of interest.
    using spline_domain_type = ddc::DiscreteDomain<bsplines_type1, bsplines_type2>;

    /**
     * @brief The type of the batch domain (obtained by removing the dimensions of interest
     * from the whole domain).
     */
    using batch_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_remove_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<evaluation_mesh_type1, evaluation_mesh_type2>>>;

    /**
     * @brief The type of the whole spline domain (cartesian product of 2D spline domain
     * and batch domain) preserving the underlying memory layout (order of dimensions).
     */
    using batched_spline_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_replace_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<evaluation_mesh_type1, evaluation_mesh_type2>,
                    ddc::detail::TypeSeq<bsplines_type1, bsplines_type2>>>;

    /// @brief The type of the extrapolation rule at the lower boundary along the first dimension.
    using left_extrapolation_rule_1_type = LeftExtrapolationRule1;

    /// @brief The type of the extrapolation rule at the upper boundary along the first dimension.
    using right_extrapolation_rule_1_type = RightExtrapolationRule1;

    /// @brief The type of the extrapolation rule at the lower boundary along the second dimension.
    using left_extrapolation_rule_2_type = LeftExtrapolationRule2;

    /// @brief The type of the extrapolation rule at the upper boundary along the second dimension.
    using right_extrapolation_rule_2_type = RightExtrapolationRule2;

private:
    LeftExtrapolationRule1 m_left_extrap_rule_1;

    RightExtrapolationRule1 m_right_extrap_rule_1;

    LeftExtrapolationRule2 m_left_extrap_rule_2;

    RightExtrapolationRule2 m_right_extrap_rule_2;

public:
    static_assert(
            std::is_same_v<LeftExtrapolationRule1,
                            typename ddc::PeriodicExtrapolationRule<
                                    tag_type1>> == bsplines_type1::is_periodic()
                    && std::is_same_v<
                               RightExtrapolationRule1,
                               typename ddc::PeriodicExtrapolationRule<
                                       tag_type1>> == bsplines_type1::is_periodic()
                    && std::is_same_v<
                               LeftExtrapolationRule2,
                               typename ddc::PeriodicExtrapolationRule<
                                       tag_type2>> == bsplines_type2::is_periodic()
                    && std::is_same_v<
                               RightExtrapolationRule2,
                               typename ddc::PeriodicExtrapolationRule<
                                       tag_type2>> == bsplines_type2::is_periodic(),
            "PeriodicExtrapolationRule has to be used if and only if dimension is periodic");
    static_assert(
            std::is_invocable_r_v<
                    double,
                    LeftExtrapolationRule1,
                    ddc::Coordinate<tag_type1>,
                    ddc::ChunkSpan<
                            double const,
                            spline_domain_type,
                            std::experimental::layout_right,
                            memory_space>>,
            "LeftExtrapolationRule1::operator() has to be callable "
            "with usual arguments.");
    static_assert(
            std::is_invocable_r_v<
                    double,
                    RightExtrapolationRule1,
                    ddc::Coordinate<tag_type1>,
                    ddc::ChunkSpan<
                            double const,
                            spline_domain_type,
                            std::experimental::layout_right,
                            memory_space>>,
            "RightExtrapolationRule1::operator() has to be callable "
            "with usual arguments.");
    static_assert(
            std::is_invocable_r_v<
                    double,
                    LeftExtrapolationRule2,
                    ddc::Coordinate<tag_type2>,
                    ddc::ChunkSpan<
                            double const,
                            spline_domain_type,
                            std::experimental::layout_right,
                            memory_space>>,
            "LeftExtrapolationRule2::operator() has to be callable "
            "with usual arguments.");
    static_assert(
            std::is_invocable_r_v<
                    double,
                    RightExtrapolationRule2,
                    ddc::Coordinate<tag_type2>,
                    ddc::ChunkSpan<
                            double const,
                            spline_domain_type,
                            std::experimental::layout_right,
                            memory_space>>,
            "RightExtrapolationRule2::operator() has to be callable "
            "with usual arguments.");

    /**
     * @brief Build a SplineEvaluator2D acting on batched_spline_domain.
     * 
     * @param left_extrap_rule1 The extrapolation rule at the lower boundary along the first dimension.
     * @param right_extrap_rule1 The extrapolation rule at the upper boundary along the first dimension.
     * @param left_extrap_rule2 The extrapolation rule at the lower boundary along the second dimension.
     * @param right_extrap_rule2 The extrapolation rule at the upper boundary along the second dimension.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    explicit SplineEvaluator2D(
            LeftExtrapolationRule1 const& left_extrap_rule1,
            RightExtrapolationRule1 const& right_extrap_rule1,
            LeftExtrapolationRule2 const& left_extrap_rule2,
            RightExtrapolationRule2 const& right_extrap_rule2)
        : m_left_extrap_rule_1(left_extrap_rule1)
        , m_right_extrap_rule_1(right_extrap_rule1)
        , m_left_extrap_rule_2(left_extrap_rule2)
        , m_right_extrap_rule_2(right_extrap_rule2)
    {
    }

    /**
     * @brief Copy-constructs.
     *
     * @param x A reference to another SplineEvaluator.
     */
    SplineEvaluator2D(SplineEvaluator2D const& x) = default;

    /**
     * @brief Move-constructs.
     *
     * @param x An rvalue to another SplineEvaluator.
     */
    SplineEvaluator2D(SplineEvaluator2D&& x) = default;

    /// @brief Destructs.
    ~SplineEvaluator2D() = default;

    /**
     * @brief Copy-assigns.
     *
     * @param x A reference to another SplineEvaluator.
     * @return A reference to this object.
     */
    SplineEvaluator2D& operator=(SplineEvaluator2D const& x) = default;

    /**
     * @brief Move-assigns.
     *
     * @param x An rvalue to another SplineEvaluator.
     * @return A reference to this object.
     */
    SplineEvaluator2D& operator=(SplineEvaluator2D&& x) = default;

    /**
     * @brief Get the lower extrapolation rule along the first dimension.
     *
     * Extrapolation rules are functors used to define the behavior of the SplineEvaluator out of the domain where the break points of the B-splines are defined.
     *
     * @return The lower extrapolation rule along the first dimension.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    left_extrapolation_rule_1_type left_extrapolation_rule_dim_1() const
    {
        return m_left_extrap_rule_1;
    }

    /**
     * @brief Get the upper extrapolation rule along the first dimension.
     *
     * Extrapolation rules are functors used to define the behavior of the SplineEvaluator out of the domain where the break points of the B-splines are defined.
     *
     * @return The upper extrapolation rule along the first dimension.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    right_extrapolation_rule_1_type right_extrapolation_rule_dim_1() const
    {
        return m_right_extrap_rule_1;
    }

    /**
     * @brief Get the lower extrapolation rule along the second dimension.
     *
     * Extrapolation rules are functors used to define the behavior of the SplineEvaluator out of the domain where the break points of the B-splines are defined.
     *
     * @return The lower extrapolation rule along the second dimension.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    left_extrapolation_rule_2_type left_extrapolation_rule_dim_2() const
    {
        return m_left_extrap_rule_2;
    }

    /**
     * @brief Get the upper extrapolation rule along the second dimension.
     *
     * Extrapolation rules are functors used to define the behavior of the SplineEvaluator out of the domain where the break points of the B-splines are defined.
     *
     * @return The upper extrapolation rule along the second dimension.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    right_extrapolation_rule_2_type right_extrapolation_rule_dim_2() const
    {
        return m_right_extrap_rule_2;
    }

    /**
     * @brief Evaluate 2D spline function (described by its spline coefficients) at a given coordinate.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * Remark: calling SplineBuilder2D then SplineEvaluator2D corresponds to a 2D spline interpolation.
     *
     * @param coord_eval The coordinate where the spline is evaluated. Note that only the components along the dimensions of interest are used.
     * @param spline_coef A ChunkSpan storing the 2D spline coefficients.
     *
     * @return The value of the spline function at the desired coordinate. 
     */
    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double operator()(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval(coord_eval, spline_coef);
    }

    /**
     * @brief Evaluate 2D spline function (described by its spline coefficients) on a mesh.
     *
     * The spline coefficients represent a 2D spline function defined on a cartesian product of batch_domain and B-splines
     * (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * This is not a nD evaluation. This is a batched 2D evaluation. This means that for each slice of coordinates
     * identified by a batch_domain_type::discrete_element_type, the evaluation is performed with the 2D set of
     * spline coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * Remark: calling SplineBuilder2D then SplineEvaluator2D corresponds to a 2D spline interpolation.
     *
     * @param[out] spline_eval The values of the 2D spline function at the desired coordinates. For practical reasons those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type.
     * @param[in] coords_eval The coordinates where the spline is evaluated. Those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 2D spline coefficients retained to perform the evaluation).
     * @param[in] spline_coef A ChunkSpan storing the 2D spline coefficients.
     */
    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void operator()(
            ddc::ChunkSpan<double, batched_evaluation_domain_type, Layout1, memory_space> const
                    spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    batched_evaluation_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, batched_spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        batch_domain_type batch_domain(coords_eval.domain());
        evaluation_domain_type1 const evaluation_domain1(spline_eval.domain());
        evaluation_domain_type2 const evaluation_domain2(spline_eval.domain());
        ddc::parallel_for_each(
                "ddc_splines_evaluate_2d",
                exec_space(),
                batch_domain,
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : evaluation_domain1) {
                        for (auto const i2 : evaluation_domain2) {
                            spline_eval_2D(i1, i2) = eval(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    /**
     * @brief Differentiate 2D spline function (described by its spline coefficients) at a given coordinate along first dimension of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be
     * obtained via various methods, such as using a SplineBuilder2D.
     *
     * @param coord_eval The coordinate where the spline is differentiated. Note that only the components along the dimensions of interest are used.
     * @param spline_coef A ChunkSpan storing the 2D spline coefficients.
     *
     * @return The derivative of the spline function at the desired coordinate. 
     */
    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv_dim_1(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval_no_bc<eval_deriv_type, eval_type>(coord_eval, spline_coef);
    }

    /**
     * @brief Differentiate 2D spline function (described by its spline coefficients) at a given coordinate along second dimension of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be
     * obtained via various methods, such as using a SplineBuilder2D.
     *
     * @param coord_eval The coordinate where the spline is differentiated. Note that only the components along the dimensions of interest are used.
     * @param spline_coef A ChunkSpan storing the 2D spline coefficients.
     *
     * @return The derivative of the spline function at the desired coordinate. 
     */
    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv_dim_2(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval_no_bc<eval_type, eval_deriv_type>(coord_eval, spline_coef);
    }

    /**
     * @brief Cross-differentiate 2D spline function (described by its spline coefficients) at a given coordinate.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be
     * obtained via various methods, such as using a SplineBuilder2D.
     *
     * @param coord_eval The coordinate where the spline is differentiated. Note that only the components along the dimensions of interest are used.
     * @param spline_coef A ChunkSpan storing the 2D spline coefficients.
     *
     * @return The derivative of the spline function at the desired coordinate. 
     */
    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv_1_and_2(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval_no_bc<eval_deriv_type, eval_deriv_type>(coord_eval, spline_coef);
    }

    /**
     * @brief Differentiate 2D spline function (described by its spline coefficients) at a given coordinate along a specified dimension of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be
     * obtained via various methods, such as using a SplineBuilder2D.
     *
     * @tparam InterestDim Dimension along which differentiation is performed.
     *
     * @param coord_eval The coordinate where the spline is differentiated. Note that only the components along the dimensions of interest are used.
     * @param spline_coef A ChunkSpan storing the 2D spline coefficients.
     *
     * @return The derivative of the spline function at the desired coordinate. 
     */
    template <class InterestDim, class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        static_assert(
                std::is_same_v<
                        InterestDim,
                        typename evaluation_mesh_type1::
                                continuous_dimension_type> || std::is_same_v<InterestDim, typename evaluation_mesh_type2::continuous_dimension_type>);
        if constexpr (std::is_same_v<
                              InterestDim,
                              typename evaluation_mesh_type1::continuous_dimension_type>) {
            return deriv_dim_1(coord_eval, spline_coef);
        } else if constexpr (std::is_same_v<
                                     InterestDim,
                                     typename evaluation_mesh_type2::continuous_dimension_type>) {
            return deriv_dim_2(coord_eval, spline_coef);
        }
    }

    /**
     * @brief Double-differentiate 2D spline function (described by its spline coefficients) at a given coordinate along specified dimensions of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be
     * obtained via various methods, such as using a SplineBuilder2D.
     *
     * Note: double-differentiation other than cross-differentiation is not supported atm. See #440
     *
     * @tparam InterestDim1 First dimension along which differentiation is performed.
     * @tparam InterestDim2 Second dimension along which differentiation is performed.
     *
     * @param coord_eval The coordinate where the spline is double-differentiated. Note that only the components along the dimensions of interest are used.
     * @param spline_coef A ChunkSpan storing the 2D spline coefficients.
     *
     * @return The derivative of the spline function at the desired coordinate. 
     */
    template <class InterestDim1, class InterestDim2, class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv2(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        static_assert(
                (std::is_same_v<
                         InterestDim1,
                         typename evaluation_mesh_type1::
                                 continuous_dimension_type> && std::is_same_v<InterestDim2, typename evaluation_mesh_type2::continuous_dimension_type>)
                || (std::is_same_v<
                            InterestDim2,
                            typename evaluation_mesh_type1::
                                    continuous_dimension_type> && std::is_same_v<InterestDim1, typename evaluation_mesh_type2::continuous_dimension_type>));
        return deriv_1_and_2(coord_eval, spline_coef);
    }

    /**
     * @brief Differentiate 2D spline function (described by its spline coefficients) on a mesh along first dimension of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a cartesian product of batch_domain and B-splines
     * (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * This is not a nD evaluation. This is a batched 2D differentiation.
     * This means that for each slice of coordinates identified by a batch_domain_type::discrete_element_type,
     * the differentiation is performed with the 2D set of spline coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * @param[out] spline_eval The derivatives of the 2D spline function at the desired coordinates. For practical reasons those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 2D spline coefficients retained to perform the evaluation).
     * @param[in] coords_eval The coordinates where the spline is differentiated. Those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 2D spline coefficients retained to perform the evaluation).
     * @param[in] spline_coef A ChunkSpan storing the 2D spline coefficients.
     */
    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv_dim_1(
            ddc::ChunkSpan<double, batched_evaluation_domain_type, Layout1, memory_space> const
                    spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    batched_evaluation_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, batched_spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        batch_domain_type batch_domain(coords_eval.domain());
        evaluation_domain_type1 const evaluation_domain1(spline_eval.domain());
        evaluation_domain_type2 const evaluation_domain2(spline_eval.domain());
        ddc::parallel_for_each(
                "ddc_splines_differentiate_2d_dim_1",
                exec_space(),
                batch_domain,
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : evaluation_domain1) {
                        for (auto const i2 : evaluation_domain2) {
                            spline_eval_2D(i1, i2) = eval_no_bc<
                                    eval_deriv_type,
                                    eval_type>(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    /**
     * @brief Differentiate 2D spline function (described by its spline coefficients) on a mesh along second dimension of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a cartesian product of batch_domain and B-splines
     * (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * This is not a nD differentiation. This is a batched 2D differentiation.
     * This means that for each slice of coordinates identified by a batch_domain_type::discrete_element_type,
     * the differentiation is performed with the 2D set of spline coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * @param[out] spline_eval The derivatives of the 2D spline function at the desired coordinates. For practical reasons those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type.
     * @param[in] coords_eval The coordinates where the spline is differentiated. Those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 2D spline coefficients retained to perform the evaluation).
     * @param[in] spline_coef A ChunkSpan storing the 2D spline coefficients.
     */
    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv_dim_2(
            ddc::ChunkSpan<double, batched_evaluation_domain_type, Layout1, memory_space> const
                    spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    batched_evaluation_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, batched_spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        batch_domain_type batch_domain(coords_eval.domain());
        evaluation_domain_type1 const evaluation_domain1(spline_eval.domain());
        evaluation_domain_type2 const evaluation_domain2(spline_eval.domain());
        ddc::parallel_for_each(
                "ddc_splines_differentiate_2d_dim_2",
                exec_space(),
                batch_domain,
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : evaluation_domain1) {
                        for (auto const i2 : evaluation_domain2) {
                            spline_eval_2D(i1, i2) = eval_no_bc<
                                    eval_type,
                                    eval_deriv_type>(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    /**
     * @brief Cross-differentiate 2D spline function (described by its spline coefficients) on a mesh along dimensions of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a cartesian product of batch_domain and B-splines
     * (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * This is not a nD cross-differentiation. This is a batched 2D cross-differentiation.
     * This means that for each slice of coordinates identified by a batch_domain_type::discrete_element_type,
     * the cross-differentiation is performed with the 2D set of spline coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * @param[out] spline_eval The cross-derivatives of the 2D spline function at the desired coordinates. For practical reasons those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type.
     * @param[in] coords_eval The coordinates where the spline is differentiated. Those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 2D spline coefficients retained to perform the evaluation).
     * @param[in] spline_coef A ChunkSpan storing the 2D spline coefficients.
     */
    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv_1_and_2(
            ddc::ChunkSpan<double, batched_evaluation_domain_type, Layout1, memory_space> const
                    spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    batched_evaluation_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, batched_spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        batch_domain_type batch_domain(coords_eval.domain());
        evaluation_domain_type1 const evaluation_domain1(spline_eval.domain());
        evaluation_domain_type2 const evaluation_domain2(spline_eval.domain());
        ddc::parallel_for_each(
                "ddc_splines_cross_differentiate",
                exec_space(),
                batch_domain,
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : evaluation_domain1) {
                        for (auto const i2 : evaluation_domain2) {
                            spline_eval_2D(i1, i2) = eval_no_bc<
                                    eval_deriv_type,
                                    eval_deriv_type>(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    /**
     * @brief Differentiate spline function (described by its spline coefficients) on a mesh along a specified dimension of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a cartesian product of batch_domain and B-splines
     * (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * This is not a nD evaluation. This is a batched 2D differentiation.
     * This means that for each slice of coordinates identified by a batch_domain_type::discrete_element_type,
     * the differentiation is performed with the 2D set of spline coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * @tparam InterestDim Dimension along which differentiation is performed.
     * @param[out] spline_eval The derivatives of the 2D spline function at the desired coordinates. For practical reasons those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type.
     * @param[in] coords_eval The coordinates where the spline is differentiated. Those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 2D spline coefficients retained to perform the evaluation).
     * @param[in] spline_coef A ChunkSpan storing the 2D spline coefficients.
     */
    template <class InterestDim, class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv(
            ddc::ChunkSpan<double, batched_evaluation_domain_type, Layout1, memory_space> const
                    spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    batched_evaluation_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, batched_spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        static_assert(
                std::is_same_v<
                        InterestDim,
                        typename evaluation_mesh_type1::
                                continuous_dimension_type> || std::is_same_v<InterestDim, typename evaluation_mesh_type2::continuous_dimension_type>);
        if constexpr (std::is_same_v<
                              InterestDim,
                              typename evaluation_mesh_type1::continuous_dimension_type>) {
            return deriv_dim_1(spline_eval, coords_eval, spline_coef);
        } else if constexpr (std::is_same_v<
                                     InterestDim,
                                     typename evaluation_mesh_type2::continuous_dimension_type>) {
            return deriv_dim_2(spline_eval, coords_eval, spline_coef);
        }
    }

    /**
     * @brief Double-differentiate 2D spline function (described by its spline coefficients) on a mesh along specified dimensions of interest.
     *
     * The spline coefficients represent a 2D spline function defined on a cartesian product of batch_domain and B-splines
     * (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * This is not a nD evaluation. This is a batched 2D differentiation.
     * This means that for each slice of coordinates identified by a batch_domain_type::discrete_element_type,
     * the differentiation is performed with the 2D set of spline coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * Note: double-differentiation other than cross-differentiation is not supported atm. See #440
     *
     * @tparam InterestDim1 First dimension along which differentiation is performed.
     * @tparam InterestDim2 Second dimension along which differentiation is performed.
     *
     * @param[out] spline_eval The derivatives of the 2D spline function at the desired coordinates. For practical reasons those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type.
     * @param[in] coords_eval The coordinates where the spline is differentiated. Those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 2D spline coefficients retained to perform the evaluation).
     * @param[in] spline_coef A ChunkSpan storing the 2D spline coefficients.
     */
    template <
            class InterestDim1,
            class InterestDim2,
            class Layout1,
            class Layout2,
            class Layout3,
            class... CoordsDims>
    void deriv2(
            ddc::ChunkSpan<double, batched_evaluation_domain_type, Layout1, memory_space> const
                    spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    batched_evaluation_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, batched_spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        static_assert(
                (std::is_same_v<
                         InterestDim1,
                         typename evaluation_mesh_type1::
                                 continuous_dimension_type> && std::is_same_v<InterestDim2, typename evaluation_mesh_type2::continuous_dimension_type>)
                || (std::is_same_v<
                            InterestDim2,
                            typename evaluation_mesh_type1::
                                    continuous_dimension_type> && std::is_same_v<InterestDim1, typename evaluation_mesh_type2::continuous_dimension_type>));
        return deriv_1_and_2(spline_eval, coords_eval, spline_coef);
    }

    /** @brief Perform batched 2D integrations of a spline function (described by its spline coefficients) along the dimensions of interest and store results on a subdomain of batch_domain.
     *
     * The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D.
     *
     * This is not a nD integration. This is a batched 2D integration.
     * This means that for each element of integrals, the integration is performed with the 2D set of
     * spline coefficients identified by the same DiscreteElement.
     *
     * @param[out] integrals The integrals of the 2D spline function on the subdomain of batch_domain. For practical reasons those are
     * stored in a ChunkSpan defined on a batch_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant.
     * @param[in] spline_coef A ChunkSpan storing the 2D spline coefficients.
     */
    template <class Layout1, class Layout2>
    void integrate(
            ddc::ChunkSpan<double, batch_domain_type, Layout1, memory_space> const integrals,
            ddc::ChunkSpan<double const, batched_spline_domain_type, Layout2, memory_space> const
                    spline_coef) const
    {
        batch_domain_type batch_domain(integrals.domain());
        ddc::Chunk values1_alloc(
                ddc::DiscreteDomain<bsplines_type1>(spline_coef.domain()),
                ddc::KokkosAllocator<double, memory_space>());
        ddc::ChunkSpan values1 = values1_alloc.span_view();
        ddc::Chunk values2_alloc(
                ddc::DiscreteDomain<bsplines_type2>(spline_coef.domain()),
                ddc::KokkosAllocator<double, memory_space>());
        ddc::ChunkSpan values2 = values2_alloc.span_view();
        Kokkos::parallel_for(
                "ddc_splines_integrate_bsplines_2d",
                Kokkos::RangePolicy<exec_space>(0, 1),
                KOKKOS_LAMBDA(int) {
                    ddc::discrete_space<bsplines_type1>().integrals(values1);
                    ddc::discrete_space<bsplines_type2>().integrals(values2);
                });

        ddc::parallel_for_each(
                "ddc_splines_integrate_bsplines",
                exec_space(),
                batch_domain,
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    integrals(j) = 0;
                    for (typename spline_domain_type1::discrete_element_type const i1 :
                         values1.domain()) {
                        for (typename spline_domain_type2::discrete_element_type const i2 :
                             values2.domain()) {
                            integrals(j) += spline_coef(i1, i2, j) * values1(i1) * values2(i2);
                        }
                    }
                });
    }

private:
    /**
     * @brief Evaluate the function on B-splines at the coordinate given.
     *
     * This function firstly deals with the boundary conditions and calls the SplineEvaluator2D::eval_no_bc function
     * to evaluate.
     *
     * @param[in] coord_eval
     * 			The 2D coordinate where we want to evaluate.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     * @param[out] vals1
     * 			A ChunkSpan with the not-null values of each function of the spline in the first dimension.
     * @param[out] vals2
     * 			A ChunkSpan with the not-null values of each function of the spline in the second dimension.
     *
     * @return A double with the value of the function at the coordinate given.
     *
     * @see SplineBoundaryValue
     */
    template <class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION double eval(
            ddc::Coordinate<CoordsDims...> coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        using Dim1 = typename evaluation_mesh_type1::continuous_dimension_type;
        using Dim2 = typename evaluation_mesh_type2::continuous_dimension_type;
        if constexpr (bsplines_type1::is_periodic()) {
            if (ddc::get<Dim1>(coord_eval) < ddc::discrete_space<bsplines_type1>().rmin()
                || ddc::get<Dim1>(coord_eval) > ddc::discrete_space<bsplines_type1>().rmax()) {
                ddc::get<Dim1>(coord_eval)
                        -= Kokkos::floor(
                                   (ddc::get<Dim1>(coord_eval)
                                    - ddc::discrete_space<bsplines_type1>().rmin())
                                   / ddc::discrete_space<bsplines_type1>().length())
                           * ddc::discrete_space<bsplines_type1>().length();
            }
        }
        if constexpr (bsplines_type2::is_periodic()) {
            if (ddc::get<Dim2>(coord_eval) < ddc::discrete_space<bsplines_type2>().rmin()
                || ddc::get<Dim2>(coord_eval) > ddc::discrete_space<bsplines_type2>().rmax()) {
                ddc::get<Dim2>(coord_eval)
                        -= Kokkos::floor(
                                   (ddc::get<Dim2>(coord_eval)
                                    - ddc::discrete_space<bsplines_type2>().rmin())
                                   / ddc::discrete_space<bsplines_type2>().length())
                           * ddc::discrete_space<bsplines_type2>().length();
            }
        }
        if constexpr (!bsplines_type1::is_periodic()) {
            if (ddc::get<Dim1>(coord_eval) < ddc::discrete_space<bsplines_type1>().rmin()) {
                return m_left_extrap_rule_1(coord_eval, spline_coef);
            }
            if (ddc::get<Dim1>(coord_eval) > ddc::discrete_space<bsplines_type1>().rmax()) {
                return m_right_extrap_rule_1(coord_eval, spline_coef);
            }
        }
        if constexpr (!bsplines_type2::is_periodic()) {
            if (ddc::get<Dim2>(coord_eval) < ddc::discrete_space<bsplines_type2>().rmin()) {
                return m_left_extrap_rule_2(coord_eval, spline_coef);
            }
            if (ddc::get<Dim2>(coord_eval) > ddc::discrete_space<bsplines_type2>().rmax()) {
                return m_right_extrap_rule_2(coord_eval, spline_coef);
            }
        }
        return eval_no_bc<eval_type, eval_type>(
                ddc::Coordinate<
                        typename evaluation_mesh_type1::continuous_dimension_type,
                        typename evaluation_mesh_type2::continuous_dimension_type>(
                        ddc::get<Dim1>(coord_eval),
                        ddc::get<Dim2>(coord_eval)),
                spline_coef);
    }

    /**
     * @brief Evaluate the function or its derivative at the coordinate given.
     *
     * @param[in] coord_eval
     * 			The coordinate where we want to evaluate.
     * @param[in] splne_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     * @tparam EvalType1
     * 			A flag indicating if we evaluate the function or its derivative in the first dimension.
     * 			The type of this object is either `eval_type` or `eval_deriv_type`.
     * @tparam EvalType2
     * 			A flag indicating if we evaluate the function or its derivative in the second dimension.
     *          The type of this object is either `eval_type` or `eval_deriv_type`.
     */
    template <class EvalType1, class EvalType2, class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION double eval_no_bc(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        static_assert(
                std::is_same_v<EvalType1, eval_type> || std::is_same_v<EvalType1, eval_deriv_type>);
        static_assert(
                std::is_same_v<EvalType2, eval_type> || std::is_same_v<EvalType2, eval_deriv_type>);
        ddc::DiscreteElement<bsplines_type1> jmin1;
        ddc::DiscreteElement<bsplines_type2> jmin2;

        std::array<double, bsplines_type1::degree() + 1> vals1_ptr;
        std::experimental::mdspan<
                double,
                std::experimental::extents<std::size_t, bsplines_type1::degree() + 1>> const
                vals1(vals1_ptr.data());
        std::array<double, bsplines_type2::degree() + 1> vals2_ptr;
        std::experimental::mdspan<
                double,
                std::experimental::extents<std::size_t, bsplines_type2::degree() + 1>> const
                vals2(vals2_ptr.data());
        ddc::Coordinate<typename evaluation_mesh_type1::continuous_dimension_type>
                coord_eval_interest1
                = ddc::select<typename evaluation_mesh_type1::continuous_dimension_type>(
                        coord_eval);
        ddc::Coordinate<typename evaluation_mesh_type2::continuous_dimension_type>
                coord_eval_interest2
                = ddc::select<typename evaluation_mesh_type2::continuous_dimension_type>(
                        coord_eval);

        if constexpr (std::is_same_v<EvalType1, eval_type>) {
            jmin1 = ddc::discrete_space<bsplines_type1>().eval_basis(vals1, coord_eval_interest1);
        } else if constexpr (std::is_same_v<EvalType1, eval_deriv_type>) {
            jmin1 = ddc::discrete_space<bsplines_type1>().eval_deriv(vals1, coord_eval_interest1);
        }
        if constexpr (std::is_same_v<EvalType2, eval_type>) {
            jmin2 = ddc::discrete_space<bsplines_type2>().eval_basis(vals2, coord_eval_interest2);
        } else if constexpr (std::is_same_v<EvalType2, eval_deriv_type>) {
            jmin2 = ddc::discrete_space<bsplines_type2>().eval_deriv(vals2, coord_eval_interest2);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_type1::degree() + 1; ++i) {
            for (std::size_t j = 0; j < bsplines_type2::degree() + 1; ++j) {
                y += spline_coef(ddc::DiscreteElement<
                                 bsplines_type1,
                                 bsplines_type2>(jmin1 + i, jmin2 + j))
                     * vals1[i] * vals2[j];
            }
        }
        return y;
    }
};
} // namespace ddc
