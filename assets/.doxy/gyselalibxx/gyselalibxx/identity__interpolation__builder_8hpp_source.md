

# File identity\_interpolation\_builder.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**identity\_interpolation\_builder.hpp**](identity__interpolation__builder_8hpp.md)

[Go to the documentation of this file](identity__interpolation__builder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_aliases.hpp"

template <class ExecSpace, class MemorySpace, class InterpolationDDim, class Basis>
class IdentityInterpolationBuilder
{
public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using continuous_dimension_type = typename InterpolationDDim::continuous_dimension_type;

    using interpolation_discrete_dimension_type = InterpolationDDim;
    using interpolation_domain_type = IdxRange<interpolation_discrete_dimension_type>;

    using basis_domain_type = typename Basis::template Impl<Basis, MemorySpace>::knot_grid;

    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_interpolation_domain_type = BatchedInterpolationGrid;

    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batch_domain_type = ddc::
            remove_dims_of_t<BatchedInterpolationGrid, interpolation_discrete_dimension_type>;

    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_basis_domain_type = ddc::replace_dim_of_t<
            BatchedInterpolationGrid,
            interpolation_discrete_dimension_type,
            basis_domain_type>;

    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_derivs_domain_type = ddc::
            remove_dims_of_t<BatchedInterpolationGrid, interpolation_discrete_dimension_type>;

public:
    static constexpr int s_nbc_xmin = 0;

    static constexpr int s_nbc_xmax = 0;

public:
    IdentityInterpolationBuilder() = default;

    template <class DataType, class Layout, class BatchedInterpolationGrid>
    void operator()(
            Field<DataType,
                  batched_basis_domain_type<BatchedInterpolationGrid>,
                  memory_space,
                  Layout> coeffs,
            ConstField<DataType, BatchedInterpolationGrid, memory_space, Layout> vals,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout>> derivs_xmax
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
            typename BatchedInterpolationGrid::discrete_vector_type nrepeat(extended_domain.size());
            BatchedInterpolationGrid repeat_domain(get_idx_range(vals).take_first(nrepeat));
            Kokkos::deep_copy(
                    coeffs[extended_domain].allocation_kokkos_view(),
                    vals[repeat_domain].allocation_kokkos_view());
        }
    }
};
```


