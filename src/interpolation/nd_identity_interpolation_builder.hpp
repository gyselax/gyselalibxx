// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

/**
 * @brief An ND builder that copies function values directly as interpolation coefficients.
 *
 * An ND extension of IdentityInterpolationBuilder. No computation is required because
 * the interpolation coefficients equal the function values at the grid nodes (e.g. for
 * an ND Lagrange interpolation on a tensor-product grid of nodes).
 *
 * @tparam ExecSpace           The Kokkos execution space.
 * @tparam MemorySpace         The Kokkos memory space.
 * @tparam DataType            The data type of field values and coefficients.
 * @tparam IdxRangeInterpolation  The ND index range for the interpolation mesh,
 *                             of the form IdxRange<Grid1, Grid2, ...>.
 * @tparam IdxRangeBasis       The ND index range for the basis types, one per
 *                             interpolation dimension, in the same order as the
 *                             grids in IdxRangeInterpolation.
 */
template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class IdxRangeInterpolation,
        class IdxRangeBasis>
class NDIdentityInterpolationBuilder;

/// The implementation of NDIdentityInterpolationBuilder. This is separate to allow a variadic Basis.
template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class IdxRangeInterpolation,
        class... Basis>
class NDIdentityInterpolationBuilder<
        ExecSpace,
        MemorySpace,
        DataType,
        IdxRangeInterpolation,
        IdxRange<Basis...>>
{
    static_assert(ddc::is_discrete_domain_v<IdxRangeInterpolation>);
    static_assert(IdxRangeInterpolation::rank() == sizeof...(Basis));
    static_assert(IdxRangeInterpolation::rank() > 0);

public:
    /// @brief The type of the Kokkos execution space.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space.
    using memory_space = MemorySpace;

    /// @brief The data type that the data is saved on.
    using data_type = DataType;

    /// @brief The ND index range for the interpolation mesh.
    using interpolation_idx_range_type = IdxRangeInterpolation;

    /// @brief The type of the index range for the bases over which coefficients of an ND Lagrange interpolation are defined.
    using coeff_idx_range_type
            = IdxRange<typename Basis::template Impl<Basis, MemorySpace>::knot_grid...>;

    /**
     * @brief Batched domain with each interpolation grid replaced by its basis domain.
     *
     * Chains ddc::replace_dim_of_t for each (InterpolationGrid_i, BasisDomain_i) pair.
     */
    template <class BatchedInterpolationIdxRange>
    using batched_basis_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                    ddc::to_type_seq_t<IdxRangeInterpolation>,
                    ddc::detail::TypeSeq<
                            typename Basis::template Impl<Basis, MemorySpace>::knot_grid...>>>;

    /// @brief The type of the index range on which derivatives should be provided (here unused).
    template <class BatchedInterpolationIdxRange>
    using batched_derivs_idx_range_type = BatchedInterpolationIdxRange;

public:
    NDIdentityInterpolationBuilder() = default;

    /**
     * @brief Compute the interpolation coefficients for a function.
     *
     * Copies vals directly to coeffs. No computation is needed because the interpolation
     * coefficients equal the function values at the grid nodes.
     *
     * @note Periodic bases are not yet supported. The copy logic for periodic dimensions
     *       (wrapping the last point back to the start) is non-trivial in ND because
     *       each combination of periodic/non-periodic dimensions must be handled
     *       separately (2^k copies for k periodic dimensions). This will be implemented
     *       when needed.
     *
     * @param[out] coeffs The coefficients of the interpolation.
     * @param[in]  vals   The values of the function on the interpolation mesh.
     */
    template <class BatchedInterpolationIdxRange>
    void operator()(
            Field<DataType,
                  batched_basis_idx_range_type<BatchedInterpolationIdxRange>,
                  memory_space> coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals) const
    {
        using IdxRangeFull = batched_basis_idx_range_type<BatchedInterpolationIdxRange>;
        using IdxRangeBatch
                = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_remove_t<
                        ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                        ddc::to_type_seq_t<IdxRangeInterpolation>>>;
        IdxRangeBatch idx_range_batch(get_idx_range(coeffs));
        coeff_idx_range_type idx_range_without_repeats(
                (ddc::discrete_space<Basis>().break_point_domain().remove_last(
                        IdxStep<typename Basis::template Impl<Basis, MemorySpace>::knot_grid>(
                                static_cast<int>(Basis::is_periodic()))))...);
        IdxRangeFull idx_range_to_fill(idx_range_batch, idx_range_without_repeats);
        Kokkos::deep_copy(
                coeffs[idx_range_to_fill].allocation_kokkos_view(),
                vals.allocation_kokkos_view());
        copy_periodic_data(coeffs, idx_range_to_fill);
    }

    /**
     * @brief Get the batched basis index range for a given batched interpolation domain.
     *
     * @param batched_interpolation_domain The full batched interpolation domain.
     * @return The batched basis index range.
     */
    template <class BatchedInterpolationIdxRange>
    batched_basis_idx_range_type<BatchedInterpolationIdxRange> batched_basis_idx_range(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        return batched_basis_idx_range_type<BatchedInterpolationIdxRange>(
                ddc::discrete_space<Basis>().full_domain()...,
                batched_interpolation_domain);
    }

private:
    template <std::size_t DimIdx = 0, class IdxRangeFull, class CoeffField>
    void copy_periodic_data(CoeffField coeffs, IdxRangeFull filled_index_range) const
    {
        using CurrentGrid = ddc::type_seq_element_t < DimIdx,
              ddc::to_type_seq_t<IdxRangeInterpolation>;
        IdxRangeFull new_filled_index_range
                = copy_periodic_data_on_dim<CurrentGrid>(coeffs, filled_index_range);
        if constexpr (DimIdx + 1 < IdxRangeInterpolation::rank()) {
            copy_periodic_data<DimIdx + 1>(coeffs, new_filled_index_range);
        }
    }

    template <class ChosenGridDim, class IdxRangeFull, class CoeffField>
    IdxRangeFull copy_periodic_data_on_dim(CoeffField coeffs, IdxRangeFull filled_index_range) const
    {
        using Dim = typename ChosenGridDim::continuous_dimension_type;
        if constexpr (Dim::PERIODIC) {
            using ChosenBasis = find_grid_t<Dim, ddc::detail::TypeSeq<Basis...>>;
            using BasisKnots =
                    typename ChosenBasis::template Impl<ChosenBasis, MemorySpace>::knot_grid;
            static_assert(ddc::in_tags_v<BasisKnots, ddc::to_type_seq_t<IdxRangeFull>>);
            using IdxRangeOthers = ddc::remove_dims_of_t<IdxRangeFull, BasisKnots>;
            IdxRangeOthers batch_idx_range(filled_index_range);
            IdxRange<BasisKnots> bp_idx_range
                    = ddc::discrete_space<ChosenBasis>().break_point_domain().remove_last(
                            IdxStep<BasisKnots>(static_cast<int>(ChosenBasis::is_periodic())));
            IdxRange<BasisKnots> extended_domain(
                    ddc::discrete_space<ChosenBasis>().full_domain().remove_first(
                            bp_idx_range.extents()));
            IdxStep<BasisKnots> nrepeat(extended_domain.size());
            IdxRange<BasisKnots> repeat_domain(
                    ddc::select<BasisKnots>(filled_index_range).take_first(nrepeat));
            IdxRangeFull to_fill(extended_domain, batch_idx_range);
            IdxRangeFull to_copy(repeat_domain, batch_idx_range);
            Kokkos::deep_copy(
                    coeffs[to_fill].allocation_kokkos_view(),
                    coeffs[to_copy].allocation_kokkos_view());
            IdxRangeFull new_filled_index_range(batch_idx_range, extended_domain);
            return new_filled_index_range;
        } else {
            return filled_index_range;
        }
    }
};
