// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cassert>
#include <type_traits>
#include <utility>

#include <experimental/mdspan>

#include <Kokkos_Core.hpp>

#include "ddc/chunk_traits.hpp"
#include "ddc/detail/macros.hpp"
#include "ddc/discrete_domain.hpp"

namespace ddc {

/** Access the domain (or subdomain) of a view
 * @param[in]  chunk the view whose domain to access
 * @return the domain of view in the queried dimensions
 */
template <class... QueryDDims, class ChunkType>
KOKKOS_FUNCTION auto get_domain(ChunkType const& chunk) noexcept
{
    static_assert(is_chunk_v<ChunkType>, "Not a chunk span type");
    return chunk.template domain<QueryDDims...>();
}

template <class ElementType, class SupportType, class LayoutStridedPolicy>
class ChunkCommon;

template <class ElementType, class... DDims, class LayoutStridedPolicy>
class ChunkCommon<ElementType, DiscreteDomain<DDims...>, LayoutStridedPolicy>
{
protected:
    /// the raw mdspan underlying this, with the same indexing (0 might no be dereferenceable)
    using internal_mdspan_type = std::experimental::mdspan<
            ElementType,
            std::experimental::dextents<std::size_t, sizeof...(DDims)>,
            std::experimental::layout_stride>;

public:
    using mdomain_type = DiscreteDomain<DDims...>;

    /// The dereferenceable part of the co-domain but with a different domain, starting at 0
    using allocation_mdspan_type = std::experimental::mdspan<
            ElementType,
            std::experimental::dextents<std::size_t, sizeof...(DDims)>,
            LayoutStridedPolicy>;

    using const_allocation_mdspan_type = std::experimental::mdspan<
            const ElementType,
            std::experimental::dextents<std::size_t, sizeof...(DDims)>,
            LayoutStridedPolicy>;

    using discrete_element_type = typename mdomain_type::discrete_element_type;

    using extents_type = typename allocation_mdspan_type::extents_type;

    using layout_type = typename allocation_mdspan_type::layout_type;

    using accessor_type = typename allocation_mdspan_type::accessor_type;

    using mapping_type = typename allocation_mdspan_type::mapping_type;

    using element_type = typename allocation_mdspan_type::element_type;

    using value_type = typename allocation_mdspan_type::value_type;

    using size_type = typename allocation_mdspan_type::size_type;

    using data_handle_type = typename allocation_mdspan_type::data_handle_type;

    using reference = typename allocation_mdspan_type::reference;

    // ChunkCommon, ChunkSpan and Chunk need to access to m_internal_mdspan and m_domain of other template versions
    template <class, class, class>
    friend class ChunkCommon;

    template <class, class, class, class>
    friend class ChunkSpan;

    template <class, class, class>
    friend class Chunk;

    static_assert(mapping_type::is_always_strided());

protected:
    /// The raw view of the data
    internal_mdspan_type m_internal_mdspan;

    /// The mesh on which this chunk is defined
    mdomain_type m_domain;

public:
    static KOKKOS_FUNCTION constexpr int rank() noexcept
    {
        return extents_type::rank();
    }

    static KOKKOS_FUNCTION constexpr int rank_dynamic() noexcept
    {
        return extents_type::rank_dynamic();
    }

    static KOKKOS_FUNCTION constexpr size_type static_extent(std::size_t r) noexcept
    {
        return extents_type::static_extent(r);
    }

    static KOKKOS_FUNCTION constexpr bool is_always_unique() noexcept
    {
        return mapping_type::is_always_unique();
    }

    static KOKKOS_FUNCTION constexpr bool is_always_exhaustive() noexcept
    {
        return mapping_type::is_always_exhaustive();
    }

    static KOKKOS_FUNCTION constexpr bool is_always_strided() noexcept
    {
        return mapping_type::is_always_strided();
    }

private:
    template <class Mapping = mapping_type>
    static KOKKOS_FUNCTION constexpr std::
            enable_if_t<std::is_constructible_v<Mapping, extents_type>, internal_mdspan_type>
            make_internal_mdspan(ElementType* ptr, mdomain_type const& domain)
    {
        if (domain.empty()) {
            return internal_mdspan_type(
                    ptr,
                    std::experimental::layout_stride::mapping<extents_type>());
        }
        extents_type extents_r(::ddc::extents<DDims>(domain).value()...);
        mapping_type mapping_r(extents_r);

        extents_type extents_s((front<DDims>(domain) + ddc::extents<DDims>(domain)).uid()...);
        std::array<std::size_t, sizeof...(DDims)> strides_s {
                mapping_r.stride(type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)...};
        std::experimental::layout_stride::mapping<extents_type> mapping_s(extents_s, strides_s);
        return internal_mdspan_type(ptr - mapping_s(front<DDims>(domain).uid()...), mapping_s);
    }

public:
    KOKKOS_FUNCTION constexpr accessor_type accessor() const
    {
        return m_internal_mdspan.accessor();
    }

    KOKKOS_FUNCTION constexpr DiscreteVector<DDims...> extents() const noexcept
    {
        return m_domain.extents();
    }

    template <class QueryDDim>
    KOKKOS_FUNCTION constexpr size_type extent() const noexcept
    {
        return m_domain.template extent<QueryDDim>();
    }

    KOKKOS_FUNCTION constexpr size_type size() const noexcept
    {
        return allocation_mdspan().size();
    }

    KOKKOS_FUNCTION constexpr mapping_type mapping() const noexcept
    {
        return allocation_mdspan().mapping();
    }

    KOKKOS_FUNCTION constexpr bool is_unique() const noexcept
    {
        return allocation_mdspan().is_unique();
    }

    KOKKOS_FUNCTION constexpr bool is_exhaustive() const noexcept
    {
        return allocation_mdspan().is_exhaustive();
    }

    KOKKOS_FUNCTION constexpr bool is_strided() const noexcept
    {
        return allocation_mdspan().is_strided();
    }

    template <class QueryDDim>
    KOKKOS_FUNCTION constexpr size_type stride() const
    {
        return m_internal_mdspan.stride(type_seq_rank_v<QueryDDim, detail::TypeSeq<DDims...>>);
    }

    /** Provide access to the domain on which this chunk is defined
     * @return the domain on which this chunk is defined
     */
    KOKKOS_FUNCTION constexpr mdomain_type domain() const noexcept
    {
        return m_domain;
    }

    /** Provide access to the domain on which this chunk is defined
     * @return the domain on which this chunk is defined
     */
    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr DiscreteDomain<QueryDDims...> domain() const noexcept
    {
        return select<QueryDDims...>(domain());
    }

protected:
    /// Empty ChunkCommon
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkCommon() = default;

    /** Constructs a new ChunkCommon from scratch
     * @param internal_mdspan
     * @param domain
     */
    KOKKOS_FUNCTION constexpr ChunkCommon(
            internal_mdspan_type internal_mdspan,
            mdomain_type const& domain) noexcept
        : m_internal_mdspan(std::move(internal_mdspan))
        , m_domain(domain)
    {
    }

    /** Constructs a new ChunkCommon from scratch
     * @param ptr the allocation pointer to the data
     * @param domain the domain that sustains the view
     */
    template <
            class Mapping = mapping_type,
            std::enable_if_t<std::is_constructible_v<Mapping, extents_type>, int> = 0>
    KOKKOS_FUNCTION constexpr ChunkCommon(ElementType* ptr, mdomain_type const& domain)
        : m_internal_mdspan(make_internal_mdspan(ptr, domain))
        , m_domain(domain)
    {
        // Handle the case where an allocation of size 0 returns a nullptr.
        assert(domain.empty() || ((ptr != nullptr) && !domain.empty()));
    }

    /** Constructs a new ChunkCommon by copy, yields a new view to the same data
     * @param other the ChunkCommon to copy
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkCommon(ChunkCommon const& other) = default;

    /** Constructs a new ChunkCommon by move
     * @param other the ChunkCommon to move
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkCommon(ChunkCommon&& other) = default;

    KOKKOS_DEFAULTED_FUNCTION ~ChunkCommon() = default;

    /** Copy-assigns a new value to this ChunkCommon, yields a new view to the same data
     * @param other the ChunkCommon to copy
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkCommon& operator=(ChunkCommon const& other) = default;

    /** Move-assigns a new value to this ChunkCommon
     * @param other the ChunkCommon to move
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr ChunkCommon& operator=(ChunkCommon&& other) = default;

    /** Access to the underlying allocation pointer
     * @return allocation pointer
     */
    KOKKOS_FUNCTION constexpr ElementType* data_handle() const
    {
        ElementType* ptr = m_internal_mdspan.data_handle();
        if (!m_domain.empty()) {
            ptr += m_internal_mdspan.mapping()(front<DDims>(m_domain).uid()...);
        }
        return ptr;
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    KOKKOS_FUNCTION constexpr internal_mdspan_type internal_mdspan() const
    {
        return m_internal_mdspan;
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    KOKKOS_FUNCTION constexpr allocation_mdspan_type allocation_mdspan() const
    {
        DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
        extents_type extents_s(::ddc::extents<DDims>(m_domain).value()...);
        if constexpr (std::is_same_v<LayoutStridedPolicy, std::experimental::layout_stride>) {
            mapping_type map(extents_s, m_internal_mdspan.mapping().strides());
            return allocation_mdspan_type(data_handle(), map);
        } else {
            mapping_type map(extents_s);
            return allocation_mdspan_type(data_handle(), map);
        }
        DDC_IF_NVCC_THEN_POP
    }
};

} // namespace ddc
