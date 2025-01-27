// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "type_seq_tools.hpp"

/**
 * @brief Copy data from a view in one layout into a span in a transposed layout.
 *
 * Layouts are described by DDC's DiscreteDomains and two layouts are considered
 * to be a transposition of one another if both domains describe data on the same
 * physical dimensions.
 *
 * @param execution_space The execution space (Host/Device) where the code will run.
 * @param transposed_field The span describing the data object which the data will be copied into.
 * @param field_to_transpose The constant span describing the data object where the original
 *                  data is found.
 *
 * @returns The transposed_field describing the data object which the data
 *          was copied into.
 */
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

/**
 * @brief If necessary transpose data into the requested dimension ordering.
 *
 * @param[in] execution_space The execution space (Host/Device) where the code will run.
 * @param[in] src The object to be transposed.
 *
 * @returns If src is already in the correct dimension ordering, return a view on src.
 *          Otherwise return a chunk with the correct dimension ordering in which the data
 *          from src has been copied.
 */
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

/**
 * @brief Create a data object in the requested dimension ordering using as allocations as possible.
 * This function does not copy data.
 *
 * @param[in] execution_space The execution space (Host/Device) where the code will run.
 * @param[in] src The object to be transposed.
 *
 * @returns If src is already in the correct dimension ordering, return a view on src.
 *          Otherwise return a chunk with the correct dimension ordering.
 */
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
