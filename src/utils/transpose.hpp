// SPDX-License-Identifier: MIT

#pragma once

/**
 * @brief Copy data from a view in one layout into a span in a transposed layout.
 *
 * Layouts are described by DDC's DiscreteDomains and two layouts are considered
 * to be a transposition of one another if both domains describe data on the same
 * physical dimensions.
 *
 * @param execution_space The execution space (Host/Device) where the code will run.
 * @param end_span The span describing the data object which the data will be copied into.
 * @param start_view The constant span describing the data object where the original
 *                  data is found.
 *
 * @returns The end_span describing the data object which the data
 *          was copied into.
 */
template <
        class ExecSpace,
        class ElementType,
        class StartDomain,
        class StartLayoutStridedPolicy,
        class MemorySpace,
        class EndDomain,
        class EndLayoutStridedPolicy>
ddc::ChunkSpan<ElementType, EndDomain, EndLayoutStridedPolicy, MemorySpace> transpose_layout(
        ExecSpace const& execution_space,
        ddc::ChunkSpan<ElementType, EndDomain, EndLayoutStridedPolicy, MemorySpace> end_span,
        ddc::ChunkView<ElementType, StartDomain, StartLayoutStridedPolicy, MemorySpace> start_view)
{
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");
    // assert that ddc::DiscreteDomain<Dims...> is a transposed discrete domain by
    // checking that it is a subset of StartDomain and that StartDomain is a subset
    // of it.
    static_assert(ddc::type_seq_contains_v<
                  ddc::to_type_seq_t<StartDomain>,
                  ddc::to_type_seq_t<EndDomain>>);
    static_assert(ddc::type_seq_contains_v<
                  ddc::to_type_seq_t<EndDomain>,
                  ddc::to_type_seq_t<StartDomain>>);

    // Check that both views have the same domain (just reordered)
    assert(start_view.domain() == end_span.domain());

    using StartIndex = typename StartDomain::discrete_element_type;
    using EndIndex = typename EndDomain::discrete_element_type;

    ddc::parallel_for_each(
            execution_space,
            start_view.domain(),
            KOKKOS_LAMBDA(StartIndex idx) { end_span(idx) = start_view(idx); });
    return end_span;
}
