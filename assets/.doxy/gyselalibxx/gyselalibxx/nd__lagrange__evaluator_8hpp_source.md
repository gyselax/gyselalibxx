

# File nd\_lagrange\_evaluator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**nd\_lagrange\_evaluator.hpp**](nd__lagrange__evaluator_8hpp.md)

[Go to the documentation of this file](nd__lagrange__evaluator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <source_location>
#include <type_traits>

#include "lagrange_evaluator.hpp"
#include "type_seq_tools.hpp"

// Forward declaration
template <class HeadEvaluator, class... Evaluators1D>
class NDLagrangeEvaluator;

template <class HeadEvaluator, class... Evaluators1D>
class NDLagrangeEvaluator
{
    static_assert(sizeof...(Evaluators1D) > 0);
    static_assert(
            (std::is_same_v<
                     typename HeadEvaluator::exec_space,
                     typename Evaluators1D::exec_space> && ...));
    static_assert(
            (std::is_same_v<
                     typename HeadEvaluator::memory_space,
                     typename Evaluators1D::memory_space> && ...));
    static_assert(
            (std::is_same_v<
                     typename HeadEvaluator::data_type,
                     typename Evaluators1D::data_type> && ...));

private:
    // Get the type of the (N-1)D LagrangeEvaluator
    using nd_minus_1_evaluator_type = std::conditional_t<
            (sizeof...(Evaluators1D) > 1),
            NDLagrangeEvaluator<Evaluators1D...>,
            ddc::type_seq_element_t<0, ddc::detail::TypeSeq<Evaluators1D...>>>;

    using HeadBasis = typename HeadEvaluator::lagrange_basis_type;
    using HeadCoeffGrid = typename HeadEvaluator::coeff_grid_type;
    using HeadContDim = typename HeadBasis::continuous_dimension_type;

    using IdxRangeHeadEval = typename HeadEvaluator::evaluation_idx_range_type;

private:
    HeadEvaluator m_head_evaluator;
    nd_minus_1_evaluator_type m_tail_evaluator;

public:
    using exec_space = typename HeadEvaluator::exec_space;

    using memory_space = typename HeadEvaluator::memory_space;

    using data_type = typename HeadEvaluator::data_type;

    using evaluation_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<type_seq_cat_t<
                    ddc::to_type_seq_t<IdxRangeHeadEval>,
                    ddc::to_type_seq_t<typename Evaluators1D::evaluation_idx_range_type>...>>;

    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batched_evaluation_idx_range_type = BatchedInterpolationIdxRange;

    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batch_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_remove_t<
                    ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                    ddc::to_type_seq_t<evaluation_idx_range_type>>>;

    using coeff_idx_range_type = IdxRange<HeadCoeffGrid, typename Evaluators1D::coeff_grid_type...>;

    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batched_coeff_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                    ddc::to_type_seq_t<evaluation_idx_range_type>,
                    ddc::to_type_seq_t<coeff_idx_range_type>>>;

public:
    explicit NDLagrangeEvaluator(HeadEvaluator const& head_eval, Evaluators1D const&... tail_evals)
        : m_head_evaluator(head_eval)
        , m_tail_evaluator(tail_evals...)
    {
    }

    NDLagrangeEvaluator(NDLagrangeEvaluator const& x) = default;

    NDLagrangeEvaluator(NDLagrangeEvaluator&& x) = default;

    ~NDLagrangeEvaluator() = default;

    NDLagrangeEvaluator& operator=(NDLagrangeEvaluator const& x) = default;

    NDLagrangeEvaluator& operator=(NDLagrangeEvaluator&& x) = default;

    template <class Layout, class... CoordsDims, class NdCoeffIdxRange>
    KOKKOS_FUNCTION data_type operator()(
            Coord<CoordsDims...> const& coord,
            ConstField<data_type, NdCoeffIdxRange, memory_space, Layout> const& coeff) const
    {
        static_assert(
                (ddc::in_tags_v<HeadContDim, ddc::detail::TypeSeq<CoordsDims...>>)&&(
                        ddc::in_tags_v<
                                typename Evaluators1D::continuous_dimension_type,
                                ddc::detail::TypeSeq<CoordsDims...>> && ...),
                "Evaluation coordinate must contain the evaluation dimensions");

        Coord<HeadContDim> coord_head(coord);
        Idx<HeadCoeffGrid> first_lagrange_knot = m_head_evaluator.find_stencil(coord_head);

        // Evaluate the Lagrange basis functions at the (adjusted) head coordinate.
        std::array<data_type, HeadBasis::degree() + 1> vals_ptr;
        Kokkos::mdspan<data_type, Kokkos::extents<std::size_t, HeadBasis::degree() + 1>> const vals(
                vals_ptr.data());
        ddc::discrete_space<HeadBasis>().eval_basis(vals, coord_head, first_lagrange_knot);

        // Tensor-product recursion: for each stencil knot, slice the coefficient
        // array along the head dimension and delegate to the (N-1)D tail evaluator.
        // When sizeof...(Evaluators1D)==1, m_tail_evaluator is the 1D LagrangeEvaluator
        // and operator() applies its configured extrapolation rules.
        data_type result = 0.;
        for (std::size_t i = 0; i < HeadBasis::degree() + 1; ++i) {
            Idx<HeadCoeffGrid> const knot_i = first_lagrange_knot + IdxStep<HeadCoeffGrid>(i);
            result += vals[i] * m_tail_evaluator(coord, coeff[knot_i]);
        }
        return result;
    }

    template <
            class Layout1,
            class Layout2,
            class Layout3,
            class IdxRangeBatched,
            class... CoordsDims>
    void operator()(
            Field<data_type, IdxRangeBatched, memory_space, Layout1> lagrange_eval,
            ConstField<Coord<CoordsDims...>, IdxRangeBatched, memory_space, Layout2> coords_eval,
            ConstField<
                    data_type,
                    batched_coeff_idx_range_type<IdxRangeBatched>,
                    memory_space,
                    Layout3> lagrange_coef) const
    {
        using IdxFull = typename IdxRangeBatched::discrete_element_type;
        using IdxBatch = typename batch_idx_range_type<IdxRangeBatched>::discrete_element_type;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                get_idx_range(lagrange_eval),
                KOKKOS_CLASS_LAMBDA(IdxFull const full_idx) {
                    IdxBatch const batch_idx(full_idx);
                    lagrange_eval(full_idx)
                            = (*this)(coords_eval(full_idx), lagrange_coef[batch_idx]);
                });
    }

    template <class Layout1, class Layout2, class IdxRangeBatched>
    void operator()(
            Field<data_type, IdxRangeBatched, memory_space, Layout1> lagrange_eval,
            ConstField<
                    data_type,
                    batched_coeff_idx_range_type<IdxRangeBatched>,
                    memory_space,
                    Layout2> lagrange_coef) const
    {
        using IdxFull = typename IdxRangeBatched::discrete_element_type;
        using IdxEval = typename evaluation_idx_range_type::discrete_element_type;
        using IdxBatch = typename batch_idx_range_type<IdxRangeBatched>::discrete_element_type;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                get_idx_range(lagrange_eval),
                KOKKOS_CLASS_LAMBDA(IdxFull const full_idx) {
                    IdxBatch const batch_idx(full_idx);
                    IdxEval const eval_idx(full_idx);
                    lagrange_eval(full_idx)
                            = (*this)(ddc::coordinate(eval_idx), lagrange_coef[batch_idx]);
                });
    }
};
```


