

# File identity\_interpolation\_builder.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**identity\_interpolation\_builder.hpp**](identity__interpolation__builder_8hpp.md)

[Go to the documentation of this file](identity__interpolation__builder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

template <class ExecSpace, class MemorySpace, class DataType, class InterpolationGrid, class Basis>
class IdentityInterpolationBuilder
{
public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using data_type = DataType;

    using continuous_dimension_type = typename InterpolationGrid::continuous_dimension_type;

    using interpolation_grid_type = InterpolationGrid;

    using interpolation_idx_range_type = IdxRange<interpolation_grid_type>;

    using basis_domain_type = typename Basis::template Impl<Basis, MemorySpace>::knot_grid;

    using deriv_type = ddc::Deriv<continuous_dimension_type>;

    template <class BatchedInterpolationIdxRange>
    using batched_basis_idx_range_type = ddc::replace_dim_of_t<
            BatchedInterpolationIdxRange,
            interpolation_grid_type,
            basis_domain_type>;

    template <class BatchedInterpolationIdxRange>
    using batched_derivs_idx_range_type = ddc::
            replace_dim_of_t<BatchedInterpolationIdxRange, interpolation_grid_type, deriv_type>;

public:
    static constexpr int s_nbc_xmin = 0;

    static constexpr int s_nbc_xmax = 0;

public:
    IdentityInterpolationBuilder() = default;

    template <class BatchedInterpolationIdxRange>
    void operator()(
            Field<DataType,
                  batched_basis_idx_range_type<BatchedInterpolationIdxRange>,
                  memory_space> coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_idx_range_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_idx_range_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmax
            = std::nullopt) const
    {
        IdxRange<basis_domain_type> bp_idx_range
                = ddc::discrete_space<Basis>().break_point_domain().remove_last(
                        IdxStep<basis_domain_type>(static_cast<int>(Basis::is_periodic())));
        Kokkos::deep_copy(
                coeffs[bp_idx_range].allocation_kokkos_view(),
                vals.allocation_kokkos_view());
        if constexpr (Basis::is_periodic()) {
            IdxRange<basis_domain_type> extended_domain(
                    ddc::discrete_space<Basis>().full_domain().remove_first(
                            bp_idx_range.extents()));
            typename BatchedInterpolationIdxRange::discrete_vector_type nrepeat(
                    extended_domain.size());
            BatchedInterpolationIdxRange repeat_domain(get_idx_range(vals).take_first(nrepeat));
            Kokkos::deep_copy(
                    coeffs[extended_domain].allocation_kokkos_view(),
                    vals[repeat_domain].allocation_kokkos_view());
        }
    }

    template <class BatchedInterpolationIdxRange>
    batched_derivs_idx_range_type<BatchedInterpolationIdxRange> batched_derivs_xmin_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> empty_deriv_range(Idx<deriv_type>(0), IdxStep<deriv_type>(0));
        return batched_derivs_idx_range_type<
                BatchedInterpolationIdxRange>(empty_deriv_range, batched_interpolation_domain);
    }

    template <class BatchedInterpolationIdxRange>
    batched_derivs_idx_range_type<BatchedInterpolationIdxRange> batched_derivs_xmax_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> empty_deriv_range(Idx<deriv_type>(0), IdxStep<deriv_type>(0));
        return batched_derivs_idx_range_type<
                BatchedInterpolationIdxRange>(empty_deriv_range, batched_interpolation_domain);
    }

    template <class BatchedInterpolationIdxRange>
    batched_basis_idx_range_type<BatchedInterpolationIdxRange> batched_basis_idx_range(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        return batched_basis_idx_range_type<BatchedInterpolationIdxRange>(
                ddc::discrete_space<Basis>().full_domain(),
                batched_interpolation_domain);
    }
};
```


