

# File transpose.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**transpose.hpp**](transpose_8hpp.md)

[Go to the documentation of this file](transpose_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "type_seq_tools.hpp"

template <
        class ExecSpace,
        class ElementType,
        class IdxRangeOut,
        class MemorySpace,
        class IdxRangeIn,
        class LayoutStridedPolicyIn,
        class LayoutStridedPolicyOut>
Field<ElementType, IdxRangeIn, MemorySpace, LayoutStridedPolicyOut> transpose_layout(
        ExecSpace const& execution_space,
        Field<ElementType, IdxRangeIn, MemorySpace, LayoutStridedPolicyOut> transposed_field,
        ConstField<ElementType, IdxRangeOut, MemorySpace, LayoutStridedPolicyIn> field_to_transpose)
{
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");
    // assert that IdxRange<Dims...> is a transposed discrete domain by
    // checking that it is a subset of IdxRangeOut and that IdxRangeOut is a subset
    // of it.
    static_assert(ddc::type_seq_contains_v<
                  ddc::to_type_seq_t<IdxRangeOut>,
                  ddc::to_type_seq_t<IdxRangeIn>>);
    static_assert(ddc::type_seq_contains_v<
                  ddc::to_type_seq_t<IdxRangeIn>,
                  ddc::to_type_seq_t<IdxRangeOut>>);

    // Check that both views have the same domain (just reordered)
    assert(get_idx_range(field_to_transpose) == get_idx_range(transposed_field));

    IdxRangeOut idx_range(field_to_transpose.domain());

    constexpr std::size_t n_dims(ddc::type_seq_size_v<ddc::to_type_seq_t<IdxRangeOut>>);

    if constexpr (n_dims < 7) {
        using ToTransposeIndex = typename IdxRangeOut::discrete_element_type;
        ddc::parallel_for_each(
                execution_space,
                idx_range,
                KOKKOS_LAMBDA(ToTransposeIndex idx) {
                    transposed_field(idx) = field_to_transpose(idx);
                });
    } else {
        using IdxRangeParallel = ddc::detail::convert_type_seq_to_discrete_domain_t<
                type_seq_range_t<ddc::to_type_seq_t<IdxRangeOut>, 0, 6>>;
        using IdxRangeSerial = ddc::detail::convert_type_seq_to_discrete_domain_t<
                type_seq_range_t<ddc::to_type_seq_t<IdxRangeOut>, 6, n_dims>>;
        using IdxParallel = typename IdxRangeParallel::discrete_element_type;
        using IdxSerial = typename IdxRangeSerial::discrete_element_type;
        IdxRangeParallel parallel_idx_range(idx_range);
        IdxRangeSerial serial_idx_range(idx_range);
        ddc::parallel_for_each(
                execution_space,
                parallel_idx_range,
                KOKKOS_LAMBDA(IdxParallel p_idx) {
                    for (IdxSerial s_idx : serial_idx_range) {
                        transposed_field(p_idx, s_idx) = field_to_transpose(p_idx, s_idx);
                    }
                });
    }
    return transposed_field;
}

namespace ddcHelper {

template <
        class TargetDomain,
        class ElementType,
        class Domain,
        class ExecSpace,
        class MemSpace,
        class FieldLayoutType>
auto create_transpose_mirror_view_and_copy(
        ExecSpace const& execution_space,
        Field<ElementType, Domain, MemSpace, FieldLayoutType> src)
{
    static_assert(
            ddc::type_seq_same_v<ddc::to_type_seq_t<Domain>, ddc::to_type_seq_t<TargetDomain>>);
    if constexpr (std::is_same_v<TargetDomain, Domain>) {
        return get_field(src);
    } else {
        TargetDomain transposed_domain(src.domain());
        using ElemType = std::remove_const_t<ElementType>;
        FieldMem<ElemType, TargetDomain, MemSpace> field_alloc(transposed_domain);
        transpose_layout(execution_space, get_field(field_alloc), get_const_field(src));
        return field_alloc;
    }
}

template <
        class TargetDomain,
        class ElementType,
        class Domain,
        class ExecSpace,
        class MemSpace,
        class FieldLayoutType>
auto create_transpose_mirror(
        ExecSpace const& execution_space,
        Field<ElementType, Domain, MemSpace, FieldLayoutType> src)
{
    static_assert(
            ddc::type_seq_same_v<ddc::to_type_seq_t<Domain>, ddc::to_type_seq_t<TargetDomain>>);
    if constexpr (std::is_same_v<TargetDomain, Domain>) {
        return get_field(src);
    } else {
        TargetDomain transposed_domain(src.domain());
        using ElemType = std::remove_const_t<ElementType>;
        FieldMem<ElemType, TargetDomain, MemSpace> field_alloc(transposed_domain);
        return field_alloc;
    }
}

} // namespace ddcHelper
```


