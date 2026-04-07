

# File nd\_identity\_interpolation\_builder.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**nd\_identity\_interpolation\_builder.hpp**](nd__identity__interpolation__builder_8hpp.md)

[Go to the documentation of this file](nd__identity__interpolation__builder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class IdxRangeInterpolation,
        class IdxRangeBasis>
class NDIdentityInterpolationBuilder;

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
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using data_type = DataType;

    using interpolation_idx_range_type = IdxRangeInterpolation;

    using coeff_idx_range_type
            = IdxRange<typename Basis::template Impl<Basis, MemorySpace>::knot_grid...>;

    template <class BatchedInterpolationIdxRange>
    using batched_basis_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                    ddc::to_type_seq_t<IdxRangeInterpolation>,
                    ddc::detail::TypeSeq<
                            typename Basis::template Impl<Basis, MemorySpace>::knot_grid...>>>;

    template <class BatchedInterpolationIdxRange>
    using batched_derivs_idx_range_type = BatchedInterpolationIdxRange;

public:
    NDIdentityInterpolationBuilder() = default;

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
        using CurrentGrid
                = ddc::type_seq_element_t<DimIdx, ddc::to_type_seq_t<IdxRangeInterpolation>>;
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
```


