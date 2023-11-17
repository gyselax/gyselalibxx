// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include <experimental/mdspan>

#include "ddc/chunk_common.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_element.hpp"

namespace ddc {

template <class, class, class>
class Chunk;

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right,
        class MemorySpace = Kokkos::DefaultHostExecutionSpace::memory_space>
class ChunkSpan;

template <class ElementType, class SupportType, class LayoutStridedPolicy, class MemorySpace>
inline constexpr bool
        enable_chunk<ChunkSpan<ElementType, SupportType, LayoutStridedPolicy, MemorySpace>> = true;

template <class ElementType, class SupportType, class LayoutStridedPolicy, class MemorySpace>
inline constexpr bool enable_borrowed_chunk<
        ChunkSpan<ElementType, SupportType, LayoutStridedPolicy, MemorySpace>> = true;

template <class ElementType, class... DDims, class LayoutStridedPolicy, class MemorySpace>
class ChunkSpan<ElementType, DiscreteDomain<DDims...>, LayoutStridedPolicy, MemorySpace>
    : public ChunkCommon<ElementType, DiscreteDomain<DDims...>, LayoutStridedPolicy>
{
protected:
    using base_type = ChunkCommon<ElementType, DiscreteDomain<DDims...>, LayoutStridedPolicy>;

    /// the raw mdspan underlying this, with the same indexing (0 might no be dereferenceable)
    using typename base_type::internal_mdspan_type;

public:
    /// type of a span of this full chunk
    using span_type
            = ChunkSpan<ElementType, DiscreteDomain<DDims...>, LayoutStridedPolicy, MemorySpace>;

    /// type of a view of this full chunk
    using view_type = ChunkSpan<
            ElementType const,
            DiscreteDomain<DDims...>,
            LayoutStridedPolicy,
            MemorySpace>;

    using mdomain_type = DiscreteDomain<DDims...>;

    using memory_space = MemorySpace;

    /// The dereferenceable part of the co-domain but with a different domain, starting at 0
    using allocation_mdspan_type = typename base_type::allocation_mdspan_type;

    using const_allocation_mdspan_type = typename base_type::const_allocation_mdspan_type;

    using discrete_element_type = typename mdomain_type::discrete_element_type;

    using extents_type = typename base_type::extents_type;

    using layout_type = typename base_type::layout_type;

    using accessor_type = typename base_type::accessor_type;

    using mapping_type = typename base_type::mapping_type;

    using element_type = typename base_type::element_type;

    using value_type = typename base_type::value_type;

    using size_type = typename base_type::size_type;

    using data_handle_type = typename base_type::data_handle_type;

    using reference = typename base_type::reference;

    template <class, class, class, class>
    friend class ChunkSpan;

protected:
    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(DiscreteElement<ODDims...> const& c) const
    {
        DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
        if constexpr (in_tags_v<QueryDDim, detail::TypeSeq<ODDims...>>) {
            return (uid<QueryDDim>(c) - front<QueryDDim>(this->m_domain).uid());
        } else {
            return std::experimental::full_extent;
        }
        DDC_IF_NVCC_THEN_POP
    }

    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(DiscreteDomain<ODDims...> const& c) const
    {
        DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
        if constexpr (in_tags_v<QueryDDim, detail::TypeSeq<ODDims...>>) {
            return std::pair<std::size_t, std::size_t>(
                    front<QueryDDim>(c) - front<QueryDDim>(this->m_domain),
                    back<QueryDDim>(c) + 1 - front<QueryDDim>(this->m_domain));
        } else {
            return std::experimental::full_extent;
        }
        DDC_IF_NVCC_THEN_POP
    }

public:
    /// Empty ChunkSpan
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkSpan() = default;

    /** Constructs a new ChunkSpan by copy, yields a new view to the same data
     * @param other the ChunkSpan to copy
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkSpan(ChunkSpan const& other) = default;

    /** Constructs a new ChunkSpan by move
     * @param other the ChunkSpan to move
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkSpan(ChunkSpan&& other) = default;

    /** Constructs a new ChunkSpan from a Chunk, yields a new view to the same data
     * @param other the Chunk to view
     */
    template <
            class OElementType,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    KOKKOS_FUNCTION constexpr ChunkSpan(
            Chunk<OElementType, mdomain_type, Allocator>& other) noexcept
        : base_type(other.m_internal_mdspan, other.m_domain)
    {
    }

    /** Constructs a new ChunkSpan from a Chunk, yields a new view to the same data
     * @param other the Chunk to view
     */
    // Disabled by SFINAE in the case of `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    KOKKOS_FUNCTION constexpr ChunkSpan(
            Chunk<OElementType, mdomain_type, Allocator> const& other) noexcept
        : base_type(other.m_internal_mdspan, other.m_domain)
    {
    }

    /** Constructs a new ChunkSpan by copy of a chunk, yields a new view to the same data
     * @param other the ChunkSpan to move
     */
    template <class OElementType>
    KOKKOS_FUNCTION constexpr ChunkSpan(
            ChunkSpan<OElementType, mdomain_type, layout_type, MemorySpace> const& other) noexcept
        : base_type(other.m_internal_mdspan, other.m_domain)
    {
    }

    /** Constructs a new ChunkSpan from scratch
     * @param ptr the allocation pointer to the data
     * @param domain the domain that sustains the view
     */
    template <
            class Mapping = mapping_type,
            std::enable_if_t<std::is_constructible_v<Mapping, extents_type>, int> = 0>
    KOKKOS_FUNCTION constexpr ChunkSpan(ElementType* const ptr, mdomain_type const& domain)
        : base_type(ptr, domain)
    {
    }

    /** Constructs a new ChunkSpan from scratch
     * @param allocation_mdspan the allocation mdspan to the data
     * @param domain the domain that sustains the view
     */
    KOKKOS_FUNCTION constexpr ChunkSpan(
            allocation_mdspan_type allocation_mdspan,
            mdomain_type const& domain)
    {
        namespace stdex = std::experimental;
        extents_type extents_s((front<DDims>(domain) + extents<DDims>(domain)).uid()...);
        std::array<std::size_t, sizeof...(DDims)> strides_s {allocation_mdspan.mapping().stride(
                type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)...};
        stdex::layout_stride::mapping<extents_type> mapping_s(extents_s, strides_s);
        this->m_internal_mdspan = internal_mdspan_type(
                allocation_mdspan.data_handle() - mapping_s(front<DDims>(domain).uid()...),
                mapping_s);
        this->m_domain = domain;
    }

    /** Constructs a new ChunkSpan from scratch
     * @param view the Kokkos view
     * @param domain the domain that sustains the view
     */
    template <class KokkosView, class = std::enable_if_t<Kokkos::is_view<KokkosView>::value>>
    KOKKOS_FUNCTION constexpr ChunkSpan(KokkosView const& view, mdomain_type const& domain) noexcept
        : ChunkSpan(
                detail::build_mdspan(view, std::make_index_sequence<sizeof...(DDims)> {}),
                domain)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION ~ChunkSpan() = default;

    /** Copy-assigns a new value to this ChunkSpan, yields a new view to the same data
     * @param other the ChunkSpan to copy
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkSpan& operator=(ChunkSpan const& other) = default;

    /** Move-assigns a new value to this ChunkSpan
     * @param other the ChunkSpan to move
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkSpan& operator=(ChunkSpan&& other) = default;

    /** Slice out some dimensions
     */
    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](
            DiscreteElement<QueryDDims...> const& slice_spec) const
    {
        auto subview = std::experimental::
                submdspan(allocation_mdspan(), get_slicer_for<DDims>(slice_spec)...);
        using detail::TypeSeq;
        using selected_meshes = type_seq_remove_t<TypeSeq<DDims...>, TypeSeq<QueryDDims...>>;
        return ChunkSpan<
                ElementType,
                decltype(select_by_type_seq<selected_meshes>(this->m_domain)),
                typename decltype(subview)::layout_type,
                memory_space>(subview, select_by_type_seq<selected_meshes>(this->m_domain));
    }

    /** Restrict to a subdomain
     */
    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](DiscreteDomain<QueryDDims...> const& odomain) const
    {
        auto subview = std::experimental::
                submdspan(allocation_mdspan(), get_slicer_for<DDims>(odomain)...);
        return ChunkSpan<
                ElementType,
                decltype(this->m_domain.restrict(odomain)),
                typename decltype(subview)::layout_type,
                memory_space>(subview, this->m_domain.restrict(odomain));
    }

    /** Element access using a list of DiscreteElement
     * @param delems discrete elements
     * @return reference to this element
     */
    template <class... DElems>
    KOKKOS_FUNCTION constexpr reference operator()(DElems const&... delems) const noexcept
    {
        static_assert(
                sizeof...(DDims) == (0 + ... + DElems::size()),
                "Invalid number of dimensions");
        static_assert((is_discrete_element_v<DElems> && ...), "Expected DiscreteElements");
        assert(((select<DDims>(take<DDims>(delems...)) >= front<DDims>(this->m_domain)) && ...));
        assert(((select<DDims>(take<DDims>(delems...)) <= back<DDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(uid<DDims>(take<DDims>(delems...))...);
    }

    /** Access to the underlying allocation pointer
     * @return allocation pointer
     */
    KOKKOS_FUNCTION constexpr ElementType* data_handle() const
    {
        return base_type::data_handle();
    }

    /** Provide a mdspan on the memory allocation
     * @return allocation mdspan
     */
    KOKKOS_FUNCTION constexpr allocation_mdspan_type allocation_mdspan() const
    {
        return base_type::allocation_mdspan();
    }

    /** Provide a mdspan on the memory allocation
     * @return allocation mdspan
     */
    KOKKOS_FUNCTION constexpr auto allocation_kokkos_view() const
    {
        auto s = this->allocation_mdspan();
        auto kokkos_layout = detail::build_kokkos_layout(
                s.extents(),
                s.mapping(),
                std::make_index_sequence<sizeof...(DDims)> {});
        return Kokkos::View<
                detail::mdspan_to_kokkos_element_t<ElementType, sizeof...(DDims)>,
                decltype(kokkos_layout),
                MemorySpace>(s.data_handle(), kokkos_layout);
    }

    KOKKOS_FUNCTION constexpr view_type span_cview() const
    {
        return view_type(*this);
    }

    KOKKOS_FUNCTION constexpr span_type span_view() const
    {
        return *this;
    }
};

template <
        class KokkosView,
        class... DDims,
        class = std::enable_if_t<Kokkos::is_view<KokkosView>::value>>
ChunkSpan(KokkosView const& view, DiscreteDomain<DDims...> domain) -> ChunkSpan<
        detail::kokkos_to_mdspan_element_t<typename KokkosView::data_type>,
        DiscreteDomain<DDims...>,
        detail::kokkos_to_mdspan_layout_t<typename KokkosView::array_layout>,
        typename KokkosView::memory_space>;

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right,
        class MemorySpace = Kokkos::HostSpace>
using ChunkView = ChunkSpan<ElementType const, SupportType, LayoutStridedPolicy, MemorySpace>;

} // namespace ddc
