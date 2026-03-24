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

/**
 * @brief Evaluates an ND Lagrange polynomial via a tensor product of 1D evaluations.
 *
 * The evaluation is recursive. Given N 1D @c LagrangeEvaluator instances (one per
 * interpolation dimension), the operator computes
 * @f[
 *   f(x_1, \ldots, x_N)
 *   = \sum_{i=0}^{d_1} L^{(1)}_i(x_1)\;
 *     \underbrace{
 *       f^{(N-1)}\!\left(x_2, \ldots, x_N;\; c[k_i, \cdot, \ldots, \cdot]\right)
 *     }_{(N-1)\text{D evaluation on slice }c[k_i,\cdot,\ldots,\cdot]}
 * @f]
 * where @f$k_i@f$ is the i-th stencil knot along the first dimension.
 *
 * The recursion terminates when only one evaluator remains: in that case the tail
 * type is the @c LagrangeEvaluator itself and the 1D evaluation (including its
 * configured extrapolation rules) is invoked directly.
 *
 * @tparam HeadEvaluator The 1D LagrangeEvaluator for the first dimension.
 * @tparam Evaluators1D  The 1D LagrangeEvaluators for the remaining dimensions.
 */
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

    // Extract the InterpolationGrid from evaluation_idx_range_type = IdxRange<Grid>.
    using IdxRangeHeadEval = typename HeadEvaluator::evaluation_idx_range_type;

private:
    HeadEvaluator m_head_evaluator;
    nd_minus_1_evaluator_type m_tail_evaluator;

public:
    /// @brief The type of the Kokkos execution space.
    using exec_space = typename HeadEvaluator::exec_space;

    /// @brief The type of the Kokkos memory space.
    using memory_space = typename HeadEvaluator::memory_space;

    /// @brief The data type.
    using data_type = typename HeadEvaluator::data_type;

    /// @brief The type of the domain for the ND evaluation mesh used by this class.
    using evaluation_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<type_seq_cat_t<
                    ddc::to_type_seq_t<IdxRangeHeadEval>,
                    ddc::to_type_seq_t<typename Evaluators1D::evaluation_idx_range_type>...>>;

    /**
     * @brief The type of the whole domain representing evaluation points.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batched_evaluation_idx_range_type = BatchedInterpolationIdxRange;

    /**
     * @brief The type of the batch domain (obtained by removing the dimensions of interest
     * from the whole domain).
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batch_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_remove_t<
                    ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                    ddc::to_type_seq_t<evaluation_idx_range_type>>>;

    /// @brief The type of the ND Lagrange domain corresponding to the dimensions of interest.
    using coeff_idx_range_type = IdxRange<HeadCoeffGrid, typename Evaluators1D::coeff_grid_type...>;

    /**
     * @brief The type of the whole Lagrange domain (cartesian product of ND Lagrange domain
     * and batch domain) preserving the order of dimensions.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batched_coeff_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                    ddc::to_type_seq_t<evaluation_idx_range_type>,
                    ddc::to_type_seq_t<coeff_idx_range_type>>>;

public:
    /**
     * @brief Construct from a sequence of 1D LagrangeEvaluators, one per dimension.
     *
     * @param head_ev   The evaluator for the first dimension.
     * @param tail_evs  The evaluators for the remaining dimensions.
     */
    explicit NDLagrangeEvaluator(HeadEvaluator const& head_ev, Evaluators1D const&... tail_evs)
        : m_head_evaluator(head_ev)
        , m_tail_evaluator(tail_evs...)
    {
    }

    /// @brief Copy-constructs.
    NDLagrangeEvaluator(NDLagrangeEvaluator const& x) = default;

    /// @brief Move-constructs.
    NDLagrangeEvaluator(NDLagrangeEvaluator&& x) = default;

    /// @brief Destructs.
    ~NDLagrangeEvaluator() = default;

    /// @brief Copy-assigns.
    NDLagrangeEvaluator& operator=(NDLagrangeEvaluator const& x) = default;

    /// @brief Move-assigns.
    NDLagrangeEvaluator& operator=(NDLagrangeEvaluator&& x) = default;

    /**
     * @brief Evaluate the ND Lagrange polynomial at a single coordinate.
     *
     * For the head dimension this method:
     * 1. Wraps the coordinate (periodic) or clips the stencil (non-periodic).
     * 2. Computes the @c degree+1 Lagrange basis values at the head coordinate.
     * 3. For each stencil knot @c k, slices @p coeff along the head dimension and
     *    calls the (N-1)D tail evaluator on the resulting slice.
     * 4. Returns the dot product of the basis values and the recursive results.
     *
     * @param coord The ND evaluation coordinate.
     * @param coeff The ND Lagrange coefficient field (no batch dimensions).
     * @return The interpolated value.
     */
    template <class Layout, class... CoordsDims, class NdCoeffIdxRange>
    KOKKOS_FUNCTION data_type operator()(
            Coord<CoordsDims...> const& coord,
            ConstField<data_type, NdCoeffIdxRange, memory_space, Layout> const coeff) const
    {
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

    /**
     * @brief Evaluate the ND Lagrange polynomial on a mesh (with explicit coordinates).
     *
     * The evaluation is parallelised over the full (batch + ND evaluation) domain.
     * For each point the single-point @c operator() is called with the corresponding
     * coordinate and batch-sliced coefficient field.
     *
     * @param[out] lagrange_eval  The interpolated values on the full batched domain.
     * @param[in]  coords_eval   The evaluation coordinates on the full batched domain.
     * @param[in]  lagrange_coef The Lagrange coefficients on the batched coefficient domain.
     */
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
        using FullIdx = typename IdxRangeBatched::discrete_element_type;
        using BatchIdx = typename batch_idx_range_type<IdxRangeBatched>::discrete_element_type;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                get_idx_range(lagrange_eval),
                KOKKOS_CLASS_LAMBDA(FullIdx const full_idx) {
                    BatchIdx const batch_idx(full_idx);
                    lagrange_eval(full_idx)
                            = (*this)(coords_eval(full_idx), lagrange_coef[batch_idx]);
                });
    }

    /**
     * @brief Evaluate the ND Lagrange polynomial on a mesh (coordinates from the grid).
     *
     * Coordinates are derived from the evaluation grid indices via @c ddc::coordinate.
     * This variant requires all evaluation grids to be genuine discrete approximations
     * of continuous dimensions.
     *
     * @param[out] lagrange_eval  The interpolated values on the full batched domain.
     * @param[in]  lagrange_coef The Lagrange coefficients on the batched coefficient domain.
     */
    template <class Layout1, class Layout2, class IdxRangeBatched>
    void operator()(
            Field<data_type, IdxRangeBatched, memory_space, Layout1> lagrange_eval,
            ConstField<
                    data_type,
                    batched_coeff_idx_range_type<IdxRangeBatched>,
                    memory_space,
                    Layout2> lagrange_coef) const
    {
        using FullIdx = typename IdxRangeBatched::discrete_element_type;
        using EvalIdx = typename evaluation_idx_range_type::discrete_element_type;
        using BatchIdx = typename batch_idx_range_type<IdxRangeBatched>::discrete_element_type;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                get_idx_range(lagrange_eval),
                KOKKOS_CLASS_LAMBDA(FullIdx const full_idx) {
                    BatchIdx const batch_idx(full_idx);
                    EvalIdx const eval_idx(full_idx);
                    lagrange_eval(full_idx)
                            = (*this)(ddc::coordinate(eval_idx), lagrange_coef[batch_idx]);
                });
    }
};
