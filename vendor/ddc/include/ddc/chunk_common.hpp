// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include <experimental/mdspan>

#include <Kokkos_Core.hpp>

#include "ddc/detail/kokkos.hpp"
#include "ddc/detail/macros.hpp"
#include "ddc/discrete_domain.hpp"

template <class T>
inline constexpr bool enable_borrowed_chunk = false;

template <class T>
inline constexpr bool enable_chunk = false;

template <class T>
inline constexpr bool is_chunk_v = enable_chunk<std::remove_const_t<std::remove_reference_t<T>>>;

template <class T>
inline constexpr bool is_borrowed_chunk_v
        = is_chunk_v<
                  T> && (std::is_lvalue_reference_v<T> || enable_borrowed_chunk<std::remove_cv_t<std::remove_reference_t<T>>>);

template <class T>
struct chunk_traits
{
    static_assert(is_chunk_v<T>);
    using value_type = std::remove_cv_t<std::remove_pointer_t<decltype(std::declval<T>().data())>>;
    using pointer_type = decltype(std::declval<T>().data());
    using reference_type = decltype(*std::declval<T>().data());
};

template <class T>
using chunk_value_t = typename chunk_traits<T>::value_type;

template <class T>
using chunk_pointer_t = typename chunk_traits<T>::pointer_type;

template <class T>
using chunk_reference_t = typename chunk_traits<T>::reference_type;

template <class T>
inline constexpr bool is_writable_chunk_v
        = !std::is_const_v<std::remove_pointer_t<chunk_pointer_t<T>>>;

/** Access the domain (or subdomain) of a view
 * @param[in]  chunk the view whose domain to access
 * @return the domain of view in the queried dimensions
 */
template <class... QueryDDims, class ChunkType>
auto get_domain(ChunkType const& chunk) noexcept
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
            std::experimental::dextents<sizeof...(DDims)>,
            std::experimental::layout_stride>;

public:
    using mdomain_type = DiscreteDomain<DDims...>;

    /// The dereferenceable part of the co-domain but with a different domain, starting at 0
    using allocation_mdspan_type = std::experimental::
            mdspan<ElementType, std::experimental::dextents<sizeof...(DDims)>, LayoutStridedPolicy>;

    using const_allocation_mdspan_type = std::experimental::mdspan<
            const ElementType,
            std::experimental::dextents<sizeof...(DDims)>,
            LayoutStridedPolicy>;

    using discrete_element_type = typename mdomain_type::discrete_element_type;

    using extents_type = typename allocation_mdspan_type::extents_type;

    using layout_type = typename allocation_mdspan_type::layout_type;

    using accessor_type = typename allocation_mdspan_type::accessor_type;

    using mapping_type = typename allocation_mdspan_type::mapping_type;

    using element_type = typename allocation_mdspan_type::element_type;

    using value_type = typename allocation_mdspan_type::value_type;

    using size_type = typename allocation_mdspan_type::size_type;

    using difference_type = typename allocation_mdspan_type::difference_type;

    using pointer = typename allocation_mdspan_type::pointer;

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
    static constexpr int rank() noexcept
    {
        return extents_type::rank();
    }

    static constexpr int rank_dynamic() noexcept
    {
        return extents_type::rank_dynamic();
    }

    static constexpr size_type static_extent(std::size_t r) noexcept
    {
        return extents_type::static_extent(r);
    }

    static constexpr bool is_always_unique() noexcept
    {
        return mapping_type::is_always_unique();
    }

    static constexpr bool is_always_contiguous() noexcept
    {
        return mapping_type::is_always_contiguous();
    }

    static constexpr bool is_always_strided() noexcept
    {
        return mapping_type::is_always_strided();
    }

public:
    constexpr accessor_type accessor() const
    {
        return m_internal_mdspan.accessor();
    }

    constexpr DiscreteVector<DDims...> extents() const noexcept
    {
        return DiscreteVector<DDims...>(
                (m_internal_mdspan.extent(type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)
                 - front<DDims>(m_domain).uid())...);
    }

    template <class QueryDDim>
    constexpr size_type extent() const noexcept
    {
        return m_internal_mdspan.extent(type_seq_rank_v<QueryDDim, detail::TypeSeq<DDims...>>)
               - front<QueryDDim>(m_domain).uid();
    }

    constexpr size_type size() const noexcept
    {
        return allocation_mdspan().size();
    }

    constexpr mapping_type mapping() const noexcept
    {
        return allocation_mdspan().mapping();
    }

    constexpr bool is_unique() const noexcept
    {
        return allocation_mdspan().is_unique();
    }

    constexpr bool is_contiguous() const noexcept
    {
        return allocation_mdspan().is_contiguous();
    }

    constexpr bool is_strided() const noexcept
    {
        return allocation_mdspan().is_strided();
    }

    template <class QueryDDim>
    constexpr size_type stride() const
    {
        return m_internal_mdspan.stride(type_seq_rank_v<QueryDDim, detail::TypeSeq<DDims...>>);
    }

    /** Provide access to the domain on which this chunk is defined
     * @return the domain on which this chunk is defined
     */
    constexpr mdomain_type domain() const noexcept
    {
        return m_domain;
    }

    /** Provide access to the domain on which this chunk is defined
     * @return the domain on which this chunk is defined
     */
    template <class... QueryDDims>
    constexpr DiscreteDomain<QueryDDims...> domain() const noexcept
    {
        return select<QueryDDims...>(domain());
    }

protected:
    /// Empty ChunkCommon
    constexpr ChunkCommon() = default;

    /** Constructs a new ChunkCommon from scratch
     * @param internal_mdspan
     * @param domain
     */
    constexpr ChunkCommon(internal_mdspan_type internal_mdspan, mdomain_type const& domain) noexcept
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
    constexpr ChunkCommon(ElementType* ptr, mdomain_type const& domain)
    {
        namespace stdex = std::experimental;
        // Handle the case where an allocation of size 0 returns a nullptr.
        assert((domain.size() == 0) || ((ptr != nullptr) && (domain.size() != 0)));

        extents_type extents_r(::extents<DDims>(domain).value()...);
        mapping_type mapping_r(extents_r);

        extents_type extents_s((front<DDims>(domain) + ::extents<DDims>(domain)).uid()...);
        std::array<std::size_t, sizeof...(DDims)> strides_s {
                mapping_r.stride(type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)...};
        stdex::layout_stride::mapping<extents_type> mapping_s(extents_s, strides_s);

        // Pointer offset to handle non-zero indexing
        ptr -= mapping_s(front<DDims>(domain).uid()...);
        m_internal_mdspan = internal_mdspan_type(ptr, mapping_s);
        m_domain = domain;
    }

    /** Constructs a new ChunkCommon by copy, yields a new view to the same data
     * @param other the ChunkCommon to copy
     */
    constexpr ChunkCommon(ChunkCommon const& other) = default;

    /** Constructs a new ChunkCommon by move
     * @param other the ChunkCommon to move
     */
    constexpr ChunkCommon(ChunkCommon&& other) = default;

    /** Copy-assigns a new value to this ChunkCommon, yields a new view to the same data
     * @param other the ChunkCommon to copy
     * @return *this
     */
    constexpr ChunkCommon& operator=(ChunkCommon const& other) = default;

    /** Move-assigns a new value to this ChunkCommon
     * @param other the ChunkCommon to move
     * @return *this
     */
    constexpr ChunkCommon& operator=(ChunkCommon&& other) = default;

    /** Access to the underlying allocation pointer
     * @return allocation pointer
     */
    constexpr ElementType* data() const
    {
        return &m_internal_mdspan(front<DDims>(m_domain).uid()...);
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    constexpr internal_mdspan_type internal_mdspan() const
    {
        return m_internal_mdspan;
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    constexpr allocation_mdspan_type allocation_mdspan() const
    {
        DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(implicit_return_from_non_void_function)
        extents_type extents_s(::extents<DDims>(m_domain).value()...);
        if constexpr (std::is_same_v<LayoutStridedPolicy, std::experimental::layout_stride>) {
            mapping_type map(extents_s, m_internal_mdspan.mapping().strides());
            return allocation_mdspan_type(data(), map);
        } else {
            mapping_type map(extents_s);
            return allocation_mdspan_type(data(), map);
        }
        DDC_IF_NVCC_THEN_POP
    }
};
