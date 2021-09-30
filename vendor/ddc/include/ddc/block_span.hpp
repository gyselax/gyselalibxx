#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include "ddc/mdomain.hpp"
#include "ddc/product_mdomain.hpp"
#include "ddc/product_mesh.hpp"

template <class, class>
class Block;

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right>
class BlockSpan;

/** Access the domain (or subdomain) of a view
 * @param[in]  view      the view whose domain to iterate
 * @return the domain of view in the queried dimensions
 */
template <class... QueryMeshes, class BlockType>
auto get_domain(BlockType const& block) noexcept
{
    return block.template domain<QueryMeshes...>();
}

template <class ElementType, class... Meshes, class LayoutStridedPolicy>
class BlockSpan<ElementType, ProductMDomain<Meshes...>, LayoutStridedPolicy>
{
public:
    using mesh_type = ProductMesh<Meshes...>;

protected:
    /// the raw mdspan underlying this, with the same indexing (0 might no be dereferenceable)
    using internal_mdspan_type = std::experimental::mdspan<
            ElementType,
            std::experimental::dextents<mesh_type::rank()>,
            std::experimental::layout_stride>;

public:
    using mdomain_type = ProductMDomain<Meshes...>;

    /// The dereferenceable part of the co-domain but with a different domain, starting at 0
    using allocation_mdspan_type = std::experimental::mdspan<
            ElementType,
            std::experimental::dextents<mesh_type::rank()>,
            LayoutStridedPolicy>;

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
    friend class BlockSpan;

    static_assert(mapping_type::is_always_strided());

protected:
    template <class QueryMesh, class... OMeshes>
    auto get_slicer_for(MCoord<OMeshes...> const& c) const
    {
        if constexpr (in_tags_v<QueryMesh, detail::TypeSeq<OMeshes...>>) {
            return get<QueryMesh>(c) - front<QueryMesh>(m_domain);
        } else {
            return std::experimental::full_extent;
        }
    }

    template <class QueryMesh, class... OMeshes>
    auto get_slicer_for(ProductMDomain<OMeshes...> const& c) const
    {
        if constexpr (in_tags_v<QueryMesh, detail::TypeSeq<OMeshes...>>) {
            return std::pair<std::size_t, std::size_t>(
                    front<QueryMesh>(c) - front<QueryMesh>(m_domain),
                    back<QueryMesh>(c) + 1 - front<QueryMesh>(m_domain));
        } else {
            return std::experimental::full_extent;
        }
    }

    /// The raw view of the data
    internal_mdspan_type m_internal_mdspan;

    /// The mesh on which this block is defined
    mdomain_type m_domain;

public:
    /** Constructs a new BlockSpan by copy, yields a new view to the same data
     * @param other the BlockSpan to copy
     */
    inline constexpr BlockSpan(BlockSpan const& other) = default;

    /** Constructs a new BlockSpan by move
     * @param other the BlockSpan to move
     */
    inline constexpr BlockSpan(BlockSpan&& other) = default;

    /** Constructs a new BlockSpan by copy of a block, yields a new view to the same data
     * @param other the BlockSpan to move
     */
    template <class OElementType>
    inline constexpr BlockSpan(Block<mdomain_type, OElementType> const& other) noexcept
        : m_internal_mdspan(other.m_internal_mdspan)
        , m_domain(other.domain())
    {
    }

    /** Constructs a new BlockSpan by copy of a block, yields a new view to the same data
     * @param other the BlockSpan to move
     */
    template <class OElementType>
    inline constexpr BlockSpan(
            BlockSpan<OElementType, mdomain_type, layout_type> const& other) noexcept
        : m_internal_mdspan(other.m_internal_mdspan)
        , m_domain(other.domain())
    {
    }

    /** Constructs a new BlockSpan from scratch
     * @param domain the domain that sustains the view
     * @param allocation_view the allocation view to the data
     */
    inline constexpr BlockSpan(allocation_mdspan_type allocation_view, mdomain_type domain)
        : m_internal_mdspan()
        , m_domain(domain)
    {
        namespace stdex = std::experimental;
        extents_type extents_s((front<Meshes>(m_domain) + ::extents<Meshes>(m_domain))...);
        std::array<std::size_t, mesh_type::rank()> strides_s {allocation_view.mapping().stride(
                type_seq_rank_v<Meshes, detail::TypeSeq<Meshes...>>)...};
        stdex::layout_stride::mapping mapping_s(extents_s, strides_s);
        m_internal_mdspan = internal_mdspan_type(
                allocation_view.data() - mapping_s(front<Meshes>(domain)...),
                mapping_s);
    }

    template <
            class Mapping = mapping_type,
            std::enable_if_t<std::is_constructible_v<Mapping, extents_type>, int> = 0>
    inline constexpr BlockSpan(ElementType* ptr, mdomain_type domain)
        : m_internal_mdspan()
        , m_domain(domain)
    {
        namespace stdex = std::experimental;
        extents_type extents_r(::extents<Meshes>(m_domain)...);
        mapping_type mapping_r(extents_r);

        extents_type extents_s((front<Meshes>(m_domain) + ::extents<Meshes>(m_domain))...);
        std::array<std::size_t, mesh_type::rank()> strides_s {
                mapping_r.stride(type_seq_rank_v<Meshes, detail::TypeSeq<Meshes...>>)...};
        stdex::layout_stride::mapping mapping_s(extents_s, strides_s);
        m_internal_mdspan
                = internal_mdspan_type(ptr - mapping_s(front<Meshes>(domain)...), mapping_s);
    }

    /** Copy-assigns a new value to this BlockSpan, yields a new view to the same data
     * @param other the BlockSpan to copy
     * @return *this
     */
    inline constexpr BlockSpan& operator=(BlockSpan const& other) = default;

    /** Move-assigns a new value to this BlockSpan
     * @param other the BlockSpan to move
     * @return *this
     */
    inline constexpr BlockSpan& operator=(BlockSpan&& other) = default;

    /** Slice out some dimensions
     * @param slices the coordinates to
     */
    template <class... QueryMeshes>
    inline constexpr auto operator[](MCoord<QueryMeshes...> const& slice_spec) const
    {
        auto subview = std::experimental::
                submdspan(allocation_mdspan(), get_slicer_for<Meshes>(slice_spec)...);
        using detail::TypeSeq;
        using selected_meshes = type_seq_remove_t<TypeSeq<Meshes...>, TypeSeq<QueryMeshes...>>;
        return ::BlockSpan(subview, select_by_type_seq<selected_meshes>(m_domain));
    }

    /** Slice out some dimensions
     * @param slices the coordinates to
     */
    template <class... QueryMeshes>
    inline constexpr auto operator[](ProductMDomain<QueryMeshes...> const& odomain) const
    {
        auto subview = std::experimental::
                submdspan(allocation_mdspan(), get_slicer_for<Meshes>(odomain)...);
        return ::BlockSpan(subview, m_domain.restrict(odomain));
    }

    template <class... OMeshes>
    inline constexpr reference operator()(
            TaggedVector<MCoordElement, OMeshes> const&... mcoords) const noexcept
    {
        assert(((mcoords >= front<OMeshes>(m_domain)) && ...));
        return m_internal_mdspan(take_first<Meshes>(mcoords...)...);
    }

    inline constexpr reference operator()(mcoord_type const& indices) const noexcept
    {
        assert(((get<Meshes>(indices) >= front<Meshes>(m_domain)) && ...));
        return m_internal_mdspan(indices.array());
    }

    template <class QueryMesh>
    inline constexpr std::size_t ibegin() const noexcept
    {
        return front<QueryMesh>(m_domain);
    }

    template <class QueryMesh>
    inline constexpr std::size_t iend() const noexcept
    {
        return back<QueryMesh>(m_domain) + 1;
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
                (m_internal_mdspan.extent(type_seq_rank_v<Meshes, detail::TypeSeq<Meshes...>>)
                 - front<Meshes>(m_domain))...);
    }

    template <class QueryMesh>
    inline constexpr size_type extent() const noexcept
    {
        return m_internal_mdspan.extent(type_seq_rank_v<QueryMesh, detail::TypeSeq<Meshes...>>)
               - front<QueryMesh>(m_domain);
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

    template <class QueryMesh>
    inline constexpr auto stride() const
    {
        return m_internal_mdspan.stride(type_seq_rank_v<QueryMesh, detail::TypeSeq<Meshes...>>);
    }

    /** Swaps this field with another
     * @param other the Block to swap with this one
     */
    inline constexpr void swap(BlockSpan& other)
    {
        BlockSpan tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    /** Provide access to the mesh on which this block is defined
     * @return the mesh on which this block is defined
     */
    inline constexpr mesh_type const& mesh() const noexcept
    {
        return m_domain.mesh();
    }

    /** Provide access to the domain on which this block is defined
     * @return the domain on which this block is defined
     */
    inline constexpr mdomain_type domain() const noexcept
    {
        return m_domain;
    }

    /** Provide access to the domain on which this block is defined
     * @return the domain on which this block is defined
     */
    template <class... QueryMeshes>
    inline constexpr ProductMDomain<QueryMeshes...> domain() const noexcept
    {
        return select<QueryMeshes...>(domain());
    }

    inline constexpr ElementType* data() const
    {
        return &m_internal_mdspan(front<Meshes>(m_domain)...);
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr internal_mdspan_type internal_view() const
    {
        return m_internal_mdspan;
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr allocation_mdspan_type allocation_mdspan() const
    {
        mapping_type m;
        extents_type extents_s(::extents<Meshes>(m_domain)...);
        if constexpr (std::is_same_v<LayoutStridedPolicy, std::experimental::layout_stride>) {
            m = mapping_type(extents_s, m_internal_mdspan.mapping().strides());
        } else {
            m = mapping_type(extents_s);
        }
        return allocation_mdspan_type(data(), m);
    }
};

template <class ElementType, class... Meshes, class Extents, class StridedLayout>
BlockSpan(
        std::experimental::mdspan<ElementType, Extents, StridedLayout> allocation_view,
        ProductMDomain<Meshes...> domain)
        -> BlockSpan<ElementType, ProductMDomain<Meshes...>, StridedLayout>;

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right>
using BlockView = BlockSpan<ElementType const, SupportType, LayoutStridedPolicy>;
