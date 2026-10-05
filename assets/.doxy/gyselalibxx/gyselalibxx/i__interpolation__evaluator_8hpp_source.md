

# File i\_interpolation\_evaluator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**i\_interpolation\_evaluator.hpp**](i__interpolation__evaluator_8hpp.md)

[Go to the documentation of this file](i__interpolation__evaluator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <type_traits>
#include <utility>

#include "ddc_aliases.hpp"

namespace detail {

template <class IdxElem>
struct ToDerivIdx;

template <class... Grids>
struct ToDerivIdx<Idx<Grids...>>
{
    using type = Idx<ddc::Deriv<typename Grids::continuous_dimension_type>...>;
};

} // namespace detail

template <class Evaluator>
struct InterpolationEvaluatorTraits
{
    using data_type = typename Evaluator::data_type;

    using evaluation_idx_range_type = typename Evaluator::evaluation_idx_range_type;

    using coord_type
            = ddc::coordinate_of_t<typename evaluation_idx_range_type::discrete_element_type>;

    using coeff_idx_range_type = typename Evaluator::coeff_idx_range_type;

    static constexpr std::size_t rank()
    {
        return evaluation_idx_range_type::rank();
    }

    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_idx_range_type<
                    BatchedInterpolationIdxRange>;

    template <class D>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_coeff_idx_range_type<D>;
};

template <
        class ExecSpace,
        class MemorySpace,
        class... BSplines,
        class... EvaluationDDim,
        class... ExtrapolationRule>
struct InterpolationEvaluatorTraits<ddc::SplineEvaluatorND<
        ExecSpace,
        MemorySpace,
        ddc::detail::TypeSeq<BSplines...>,
        ddc::detail::TypeSeq<EvaluationDDim...>,
        ddc::detail::TypeSeq<ExtrapolationRule...>>>
{
private:
    using Evaluator = ddc::SplineEvaluatorND<
            ExecSpace,
            MemorySpace,
            ddc::detail::TypeSeq<BSplines...>,
            ddc::detail::TypeSeq<EvaluationDDim...>,
            ddc::detail::TypeSeq<ExtrapolationRule...>>;

    template <std::size_t... Is>
    static auto make_evaluation_idx_range_type(std::index_sequence<Is...>) ->
            typename Evaluator::template evaluation_domain_type<Is...>;

    template <std::size_t... Is>
    static auto make_coeff_idx_range_type(std::index_sequence<Is...>) ->
            typename Evaluator::template spline_domain_type<Is...>;

public:
    using data_type = double;

    using evaluation_idx_range_type = decltype(make_evaluation_idx_range_type(
            std::make_index_sequence<sizeof...(BSplines)>()));

    using coord_type = Coord<typename EvaluationDDim::continuous_dimension_type...>;

    using coeff_idx_range_type
            = decltype(make_coeff_idx_range_type(std::make_index_sequence<sizeof...(BSplines)>()));

    static constexpr std::size_t rank()
    {
        return sizeof...(BSplines);
    }

    template <class BatchedInterpolationIdxRange>
    using batched_evaluation_idx_range_type =
            typename Evaluator::template batched_evaluation_domain_type<
                    BatchedInterpolationIdxRange>;

    template <class BatchedInterpolationIdxRange>
    using batched_coeff_idx_range_type =
            typename Evaluator::template batched_spline_domain_type<BatchedInterpolationIdxRange>;
};


namespace concepts {

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
```


