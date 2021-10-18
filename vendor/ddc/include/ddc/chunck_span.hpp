#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include <experimental/mdspan>

#include "ddc/detail/discrete_space.hpp"
#include "ddc/discrete_domain.hpp"

template <class, class>
class Chunk;

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right>
class ChunkSpan;

/** Access the domain (or subdomain) of a view
 * @param[in]  view      the view whose domain to iterate
 * @return the domain of view in the queried dimensions
 */
template <class... QueryDDims, class ChunkType>
auto get_domain(ChunkType const& chunck) noexcept
{
    return chunck.template domain<QueryDDims...>();
}

template <class ElementType, class... DDims, class LayoutStridedPolicy>
class ChunkSpan<ElementType, DiscreteDomain<DDims...>, LayoutStridedPolicy>
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

    using mcoord_type = typename mdomain_type::mcoord_type;

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

    template <class, class, class>
    friend class ChunkSpan;

    static_assert(mapping_type::is_always_strided());

protected:
    template <class QueryDDim, class... ODDims>
    auto get_slicer_for(DiscreteCoordinate<ODDims...> const& c) const
    {
        if constexpr (in_tags_v<QueryDDim, detail::TypeSeq<ODDims...>>) {
            return get<QueryDDim>(c) - front<QueryDDim>(m_domain);
        } else {
            return std::experimental::full_extent;
        }
    }

    template <class QueryDDim, class... ODDims>
    auto get_slicer_for(DiscreteDomain<ODDims...> const& c) const
    {
        if constexpr (in_tags_v<QueryDDim, detail::TypeSeq<ODDims...>>) {
            return std::pair<std::size_t, std::size_t>(
                    front<QueryDDim>(c) - front<QueryDDim>(m_domain),
                    back<QueryDDim>(c) + 1 - front<QueryDDim>(m_domain));
        } else {
            return std::experimental::full_extent;
        }
    }

    /// The raw view of the data
    internal_mdspan_type m_internal_mdspan;

    /// The mesh on which this chunck is defined
    mdomain_type m_domain;

public:
    /** Constructs a new ChunkSpan by copy, yields a new view to the same data
     * @param other the ChunkSpan to copy
     */
    inline constexpr ChunkSpan(ChunkSpan const& other) = default;

    /** Constructs a new ChunkSpan by move
     * @param other the ChunkSpan to move
     */
    inline constexpr ChunkSpan(ChunkSpan&& other) = default;

    /** Constructs a new ChunkSpan by copy of a chunck, yields a new view to the same data
     * @param other the ChunkSpan to move
     */
    template <class OElementType>
    inline constexpr ChunkSpan(Chunk<mdomain_type, OElementType> const& other) noexcept
        : m_internal_mdspan(other.m_internal_mdspan)
        , m_domain(other.domain())
    {
    }

    /** Constructs a new ChunkSpan by copy of a chunck, yields a new view to the same data
     * @param other the ChunkSpan to move
     */
    template <class OElementType>
    inline constexpr ChunkSpan(
            ChunkSpan<OElementType, mdomain_type, layout_type> const& other) noexcept
        : m_internal_mdspan(other.m_internal_mdspan)
        , m_domain(other.domain())
    {
    }

    /** Constructs a new ChunkSpan from scratch
     * @param domain the domain that sustains the view
     * @param allocation_view the allocation view to the data
     */
    inline constexpr ChunkSpan(allocation_mdspan_type allocation_view, mdomain_type domain)
        : m_internal_mdspan()
        , m_domain(domain)
    {
        namespace stdex = std::experimental;
        extents_type extents_s((front<DDims>(m_domain) + ::extents<DDims>(m_domain))...);
        std::array<std::size_t, sizeof...(DDims)> strides_s {allocation_view.mapping().stride(
                type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)...};
        stdex::layout_stride::mapping mapping_s(extents_s, strides_s);
        m_internal_mdspan = internal_mdspan_type(
                allocation_view.data() - mapping_s(front<DDims>(domain)...),
                mapping_s);
    }

    template <
            class Mapping = mapping_type,
            std::enable_if_t<std::is_constructible_v<Mapping, extents_type>, int> = 0>
    inline constexpr ChunkSpan(ElementType* ptr, mdomain_type domain)
        : m_internal_mdspan()
        , m_domain(domain)
    {
        namespace stdex = std::experimental;
        extents_type extents_r(::extents<DDims>(m_domain)...);
        mapping_type mapping_r(extents_r);

        extents_type extents_s((front<DDims>(m_domain) + ::extents<DDims>(m_domain))...);
        std::array<std::size_t, sizeof...(DDims)> strides_s {
                mapping_r.stride(type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)...};
        stdex::layout_stride::mapping mapping_s(extents_s, strides_s);
        m_internal_mdspan
                = internal_mdspan_type(ptr - mapping_s(front<DDims>(domain)...), mapping_s);
    }

    /** Copy-assigns a new value to this ChunkSpan, yields a new view to the same data
     * @param other the ChunkSpan to copy
     * @return *this
     */
    inline constexpr ChunkSpan& operator=(ChunkSpan const& other) = default;

    /** Move-assigns a new value to this ChunkSpan
     * @param other the ChunkSpan to move
     * @return *this
     */
    inline constexpr ChunkSpan& operator=(ChunkSpan&& other) = default;

    /** Slice out some dimensions
     */
    template <class... QueryDDims>
    inline constexpr auto operator[](DiscreteCoordinate<QueryDDims...> const& slice_spec) const
    {
        auto subview = std::experimental::
                submdspan(allocation_mdspan(), get_slicer_for<DDims>(slice_spec)...);
        using detail::TypeSeq;
        using selected_meshes = type_seq_remove_t<TypeSeq<DDims...>, TypeSeq<QueryDDims...>>;
        return ::ChunkSpan(subview, select_by_type_seq<selected_meshes>(m_domain));
    }

    /** Slice out some dimensions
     */
    template <class... QueryDDims>
    inline constexpr auto operator[](DiscreteDomain<QueryDDims...> const& odomain) const
    {
        auto subview = std::experimental::
                submdspan(allocation_mdspan(), get_slicer_for<DDims>(odomain)...);
        return ::ChunkSpan(subview, m_domain.restrict(odomain));
    }

    // Warning: Do not use DiscreteCoordinate because of template deduction issue with clang 12
    template <class... ODDims>
    inline constexpr reference operator()(
            detail::TaggedVector<DiscreteCoordElement, ODDims> const&... mcoords) const noexcept
    {
        assert(((mcoords >= front<ODDims>(m_domain)) && ...));
        return m_internal_mdspan(take_first<DDims>(mcoords...)...);
    }

    inline constexpr reference operator()(mcoord_type const& indices) const noexcept
    {
        assert(((get<DDims>(indices) >= front<DDims>(m_domain)) && ...));
        return m_internal_mdspan(indices.array());
    }

    /// @deprecated
    template <class QueryDDim>
    [[deprecated]] inline constexpr std::size_t ibegin() const noexcept
    {
        return front<QueryDDim>(m_domain);
    }

    /// @deprecated
    template <class QueryDDim>
    [[deprecated]] inline constexpr std::size_t iend() const noexcept
    {
        return back<QueryDDim>(m_domain) + 1;
    }

    inline accessor_type accessor() const
    {
        return m_internal_mdspan.accessor();
    }

    static inline constexpr int rank() noexcept
    {
        return extents_type::rank();
    }

    static inline constexpr int rank_dynamic() noexcept
    {
        return extents_type::rank_dynamic();
    }

    // static inline constexpr size_type static_extent(size_t r) noexcept
    // {
    //     return extents_type::static_extent(r);
    // }

    inline constexpr mcoord_type extents() const noexcept
    {
        return mcoord_type(
                (m_internal_mdspan.extent(type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)
                 - front<DDims>(m_domain))...);
    }

    template <class QueryDDim>
    inline constexpr size_type extent() const noexcept
    {
        return m_internal_mdspan.extent(type_seq_rank_v<QueryDDim, detail::TypeSeq<DDims...>>)
               - front<QueryDDim>(m_domain);
    }

    inline constexpr size_type size() const noexcept
    {
        return allocation_mdspan().size();
    }

    inline constexpr size_type unique_size() const noexcept
    {
        return allocation_mdspan().unique_size();
    }

    static inline constexpr bool is_always_unique() noexcept
    {
        return mapping_type::is_always_unique();
    }

    static inline constexpr bool is_always_contiguous() noexcept
    {
        return mapping_type::is_always_contiguous();
    }

    static inline constexpr bool is_always_strided() noexcept
    {
        return mapping_type::is_always_strided();
    }

    inline constexpr mapping_type mapping() const noexcept
    {
        return allocation_mdspan().mapping();
    }

    inline constexpr bool is_unique() const noexcept
    {
        return allocation_mdspan().is_unique();
    }

    inline constexpr bool is_contiguous() const noexcept
    {
        return allocation_mdspan().is_contiguous();
    }

    inline constexpr bool is_strided() const noexcept
    {
        return allocation_mdspan().is_strided();
    }

    template <class QueryDDim>
    inline constexpr auto stride() const
    {
        return m_internal_mdspan.stride(type_seq_rank_v<QueryDDim, detail::TypeSeq<DDims...>>);
    }

    /** Swaps this field with another
     * @param other the Chunk to swap with this one
     */
    inline constexpr void swap(ChunkSpan& other)
    {
        ChunkSpan tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    /** Provide access to the domain on which this chunck is defined
     * @return the domain on which this chunck is defined
     */
    inline constexpr mdomain_type domain() const noexcept
    {
        return m_domain;
    }

    /** Provide access to the domain on which this chunck is defined
     * @return the domain on which this chunck is defined
     */
    template <class... QueryDDims>
    inline constexpr DiscreteDomain<QueryDDims...> domain() const noexcept
    {
        return select<QueryDDims...>(domain());
    }

    inline constexpr ElementType* data() const
    {
        return &m_internal_mdspan(front<DDims>(m_domain)...);
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr internal_mdspan_type internal_mdspan() const
    {
        return m_internal_mdspan;
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr allocation_mdspan_type allocation_mdspan() const
    {
        mapping_type m;
        extents_type extents_s(::extents<DDims>(m_domain)...);
        if constexpr (std::is_same_v<LayoutStridedPolicy, std::experimental::layout_stride>) {
            // Temporary workaround: layout_stride::mapping is missing the function `strides`
            const std::array<std::size_t, extents_type::rank()> strides {
                    m_internal_mdspan.mapping().stride(
                            type_seq_rank_v<DDims, detail::TypeSeq<DDims...>>)...};
            m = mapping_type(extents_s, strides);
        } else {
            m = mapping_type(extents_s);
        }
        return allocation_mdspan_type(data(), m);
    }
};

template <class ElementType, class... DDims, class Extents, class StridedLayout>
ChunkSpan(
        std::experimental::mdspan<ElementType, Extents, StridedLayout> allocation_view,
        DiscreteDomain<DDims...> domain)
        -> ChunkSpan<ElementType, DiscreteDomain<DDims...>, StridedLayout>;

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right>
using ChunkView = ChunkSpan<ElementType const, SupportType, LayoutStridedPolicy>;
