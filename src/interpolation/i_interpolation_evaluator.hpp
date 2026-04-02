// SPDX-License-Identifier: MIT
#pragma once

#include <type_traits>

#include "ddc_aliases.hpp"

namespace detail {

/**
 * @brief Convert an index type over grids to the matching index type over derivative dimensions.
 *
 * Given @c Idx<Grid1, Grid2, ...>, produces
 * @c Idx<ddc::Deriv<Grid1::continuous_dimension_type>, ddc::Deriv<Grid2::continuous_dimension_type>, ...>.
 *
 * @tparam IdxElem An @c Idx<Grids...> type.
 */
template <class IdxElem>
struct ToDerivIdx;

template <class... Grids>
struct ToDerivIdx<Idx<Grids...>>
{
    using type = Idx<ddc::Deriv<typename Grids::continuous_dimension_type>...>;
};

} // namespace detail

/**
 * @brief A traits struct for accessing type aliases of an interpolation evaluator.
 *
 * The primary template delegates to the evaluator's own type aliases. It handles
 * the common case where an evaluator already uses the convention names (e.g.
 * LagrangeEvaluator exposes @c evaluation_idx_range_type and
 * @c batched_coeff_idx_range_type directly).
 *
 * Specialise this struct to adapt external evaluators whose alias names differ
 * (e.g. ddc::SplineEvaluator).
 *
 * Defines:
 *   Type aliases:
 *   - data_type
 *   - evaluation_idx_range_type
 *   - coord_type
 *   - coeff_idx_range_type
 *   Static functions
 *   - rank()
 *   Type calculators
 *   - batched_evaluation_idx_range_type
 *   - batched_coeff_idx_range_type
 *
 * @tparam Evaluator The interpolation evaluator type.
 */
template <class Evaluator>
struct InterpolationEvaluatorTraits
{
    /// @brief The data type that the data is saved on.
    using data_type = typename Evaluator::data_type;

    /// @brief The ND index range for the evaluation mesh.
    using evaluation_idx_range_type = typename Evaluator::evaluation_idx_range_type;

    /// @brief The ND coordinate type corresponding to the evaluation mesh.
    using coord_type
            = ddc::coordinate_of_t<typename evaluation_idx_range_type::discrete_element_type>;

    /// @brief The type of the ND index range on which the interpolation coefficients are defined.
    using coeff_idx_range_type = typename Evaluator::coeff_idx_range_type;

    /// @brief The number of interpolation dimensions.
    static constexpr std::size_t rank()
    {
        return evaluation_idx_range_type::rank();
    }

    /// @brief Batched index range for the evaluation
    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_idx_range_type<
                    BatchedInterpolationIdxRange>;

    /// @brief Batched domain with the evaluation grid replaced by the coefficient grid(s).
    template <class D>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_coeff_idx_range_type<D>;
};

/**
 * @brief Specialisation of InterpolationEvaluatorTraits for ddc::SplineEvaluator.
 *
 * ddc::SplineEvaluator uses different alias names from the InterpolationEvaluator
 * convention. This specialisation provides the mapping so that ddc::SplineEvaluator
 * can be used directly as an InterpolationEvaluator.
 *
 * Mapping:
 *   evaluation_discrete_dimension_type -> (defines evaluation_idx_range_type)
 *   evaluation_domain_type             -> evaluation_idx_range_type
 *   spline_domain_type                 -> coeff_idx_range_type
 *   batched_spline_domain_type<D>      -> batched_coeff_idx_range_type<D>
 */
template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class EvaluationDDim,
        class LowerExtrapolationRule,
        class UpperExtrapolationRule>
struct InterpolationEvaluatorTraits<ddc::SplineEvaluator<
        ExecSpace,
        MemorySpace,
        BSplines,
        EvaluationDDim,
        LowerExtrapolationRule,
        UpperExtrapolationRule>>
{
private:
    using Evaluator = ddc::SplineEvaluator<
            ExecSpace,
            MemorySpace,
            BSplines,
            EvaluationDDim,
            LowerExtrapolationRule,
            UpperExtrapolationRule>;

public:
    /// @brief The data type that the data is saved on.
    using data_type = double;

    /// @brief The 1D index range for the evaluation mesh.
    using evaluation_idx_range_type = typename Evaluator::evaluation_domain_type;

    /// @brief The 1D coordinate type corresponding to the evaluation mesh.
    using coord_type = Coord<typename EvaluationDDim::continuous_dimension_type>;

    /// @brief The type of the ND index range on which the interpolation coefficients are defined.
    using coeff_idx_range_type = typename Evaluator::spline_domain_type;

    /// @brief The number of interpolation dimensions (always 1 for SplineEvaluator).
    static constexpr std::size_t rank()
    {
        return 1;
    }

    /// @brief Batched index range for the evaluation
    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_domain_type<
                    BatchedInterpolationIdxRange>;

    /// @brief Batched domain with the evaluation grid replaced by BSplines.
    template <class BatchedInterpolationIdxRange>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_spline_domain_type<BatchedInterpolationIdxRange>;
};

namespace concepts {

/**
 * @brief A concept describing an ND interpolation evaluator.
 *
 * An interpolation evaluator is a callable that takes a field of interpolation
 * coefficients (on a basis domain) and produces function values on an evaluation
 * mesh, optionally at user-supplied coordinates. The evaluator may operate over
 * one or more interpolation dimensions simultaneously (N ≥ 1).
 *
 * Type information is accessed through InterpolationEvaluatorTraits<Evaluator>, which
 * has a primary template delegating to the evaluator's own aliases and can be
 * specialised for external evaluators (e.g. ddc::SplineEvaluator).
 *
 * The template alias and method requirements are verified using evaluation_idx_range_type
 * as a representative domain.
 */
template <class Evaluator>
concept InterpolationEvaluator = requires
{
    typename Evaluator::exec_space;
    typename Evaluator::memory_space;
    typename InterpolationEvaluatorTraits<Evaluator>::data_type;
    typename InterpolationEvaluatorTraits<Evaluator>::evaluation_idx_range_type;
    typename InterpolationEvaluatorTraits<Evaluator>::coord_type;
    typename InterpolationEvaluatorTraits<Evaluator>::coeff_idx_range_type;
}
&&requires()
{
    {
        InterpolationEvaluatorTraits<Evaluator>::rank()
        } -> std::same_as<std::size_t>;
}
&&requires(
        Evaluator const& e,
        Field<typename InterpolationEvaluatorTraits<Evaluator>::data_type,
              typename InterpolationEvaluatorTraits<Evaluator>::evaluation_idx_range_type,
              typename Evaluator::memory_space> eval,
        ConstField<
                typename InterpolationEvaluatorTraits<Evaluator>::data_type,
                typename InterpolationEvaluatorTraits<Evaluator>::
                        template batched_coeff_idx_range_type<typename InterpolationEvaluatorTraits<
                                Evaluator>::evaluation_idx_range_type>,
                typename Evaluator::memory_space> coeffs,
        ConstField<
                typename InterpolationEvaluatorTraits<Evaluator>::coord_type,
                typename InterpolationEvaluatorTraits<Evaluator>::evaluation_idx_range_type,
                typename Evaluator::memory_space> coords,
        typename detail::ToDerivIdx<typename InterpolationEvaluatorTraits<
                Evaluator>::evaluation_idx_range_type::discrete_element_type>::type deriv_order)
{
    {e(eval, coeffs)};
    {e(eval, coords, coeffs)};
    {e.deriv(deriv_order, eval, coeffs)};
    {e.deriv(deriv_order, eval, coords, coeffs)};
}
&&requires(
        Evaluator const& e,
        ConstField<
                typename InterpolationEvaluatorTraits<Evaluator>::data_type,
                typename InterpolationEvaluatorTraits<Evaluator>::
                        template batched_coeff_idx_range_type<typename InterpolationEvaluatorTraits<
                                Evaluator>::evaluation_idx_range_type>,
                typename Evaluator::memory_space> coeffs,
        typename InterpolationEvaluatorTraits<Evaluator>::coord_type coord,
        typename detail::ToDerivIdx<typename InterpolationEvaluatorTraits<
                Evaluator>::evaluation_idx_range_type::discrete_element_type>::type deriv_order)
{
    {
        e(coord, coeffs)
        } -> std::same_as<typename InterpolationEvaluatorTraits<Evaluator>::data_type>;
    {
        e.deriv(deriv_order, coord, coeffs)
        } -> std::same_as<typename InterpolationEvaluatorTraits<Evaluator>::data_type>;
};

} // namespace concepts
