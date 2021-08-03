#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include "mdomain.h"
#include "product_mdomain.h"
#include "product_mesh.h"
#include "view.h"

template <class, class>
class Block;

template <
        class SupportType,
        class ElementType,
        class LayoutPolicy = std::experimental::layout_stride>
class BlockView;

/** Access the domain (or subdomain) of a view
 * @param[in]  view      the view whose domain to iterate
 * @return the domain of view in the queried dimensions
 */
template <class... QueryMeshes, class BlockType>
auto get_domain(BlockType const& block) noexcept
{
    return block.template domain<QueryMeshes...>();
}

template <class... Meshes, class ElementType>
class BlockView<ProductMDomain<Meshes...>, ElementType, std::experimental::layout_stride>
{
    using underlying_raw_view_type = std::experimental::mdspan<
            ElementType,
            std::experimental::dextents<sizeof...(Meshes)>,
            std::experimental::layout_right>;

public:
    using mdomain_type = ProductMDomain<Meshes...>;

    using mesh_type = ProductMesh<Meshes...>;

    /// ND memory view
    using raw_view_type = std::experimental::mdspan<
            ElementType,
            std::experimental::dextents<mesh_type::rank()>,
            std::experimental::layout_stride>;

    using mcoord_type = typename mdomain_type::mcoord_type;

    using extents_type = typename raw_view_type::extents_type;

    using layout_type = typename raw_view_type::layout_type;

    using accessor_type = typename raw_view_type::accessor_type;

    using mapping_type = typename raw_view_type::mapping_type;

    using element_type = typename raw_view_type::element_type;

    using value_type = typename raw_view_type::value_type;

    using size_type = typename raw_view_type::size_type;

    using difference_type = typename raw_view_type::difference_type;

    using pointer = typename raw_view_type::pointer;

    using reference = typename raw_view_type::reference;

    template <class, class, class>
    friend class BlockView;

protected:
    template <class QueryMesh, class... OMeshes>
    static auto get_slicer_for(MCoord<OMeshes...> const& c)
    {
        if constexpr (in_tags_v<QueryMesh, detail::TypeSeq<OMeshes...>>) {
            return get<QueryMesh>(c);
        } else {
            return std::experimental::full_extent;
        }
    }

    template <class QueryMesh, class... OMeshes>
    static auto get_slicer_for(ProductMDomain<OMeshes...> const& c)
    {
        if constexpr (in_tags_v<QueryMesh, detail::TypeSeq<OMeshes...>>) {
            return std::pair<std::size_t, std::size_t>(0, back<QueryMesh>(c) + 1);
        } else {
            return std::experimental::full_extent;
        }
    }

    /// The raw view of the data
    raw_view_type m_raw;

    /// The mesh on which this block is defined
    mdomain_type m_domain;

public:
    /** Constructs a new BlockView by copy, yields a new view to the same data
     * @param other the BlockView to copy
     */
    inline constexpr BlockView(BlockView const& other) = default;

    /** Constructs a new BlockView by move
     * @param other the BlockView to move
     */
    inline constexpr BlockView(BlockView&& other) = default;

    /** Constructs a new BlockView by copy of a block, yields a new view to the same data
     * @param other the BlockView to move
     */
    template <class OElementType>
    inline constexpr BlockView(Block<mdomain_type, OElementType> const& other) noexcept
        : m_raw(other.m_raw)
        , m_domain(other.domain())
    {
    }

    /** Constructs a new BlockView by copy of a block, yields a new view to the same data
     * @param other the BlockView to move
     */
    template <class OElementType>
    inline constexpr BlockView(
            BlockView<mdomain_type, OElementType, std::experimental::layout_stride> const&
                    other) noexcept
        : m_raw(other.m_raw)
        , m_domain(other.domain())
    {
    }

    /** Constructs a new BlockView from scratch
     * @param domain the domain that sustains the view
     * @param raw_view the raw view to the data
     */
    inline constexpr BlockView(mdomain_type domain, raw_view_type raw_view)
        : m_raw(raw_view)
        , m_domain(domain)
    {
    }

    inline constexpr BlockView(mdomain_type domain, ElementType* ptr) : m_raw(), m_domain(domain)
    {
        namespace stdex = std::experimental;
        stdex::dextents<mesh_type::rank()> extents_r(::extents<Meshes>(m_domain)...);
        stdex::layout_right::mapping mapping_r(extents_r);

        stdex::dextents<mesh_type::rank()> extents_s(
                (front<Meshes>(m_domain) + ::extents<Meshes>(m_domain))...);
        std::array<std::size_t, mesh_type::rank()> strides_s {
                mapping_r.stride(type_seq_rank_v<Meshes, detail::TypeSeq<Meshes...>>)...};
        stdex::layout_stride::mapping mapping_s(extents_s, strides_s);
        m_raw = raw_view_type(ptr - mapping_s(front<Meshes>(domain)...), mapping_s);
    }

    /** Copy-assigns a new value to this BlockView, yields a new view to the same data
     * @param other the BlockView to copy
     * @return *this
     */
    inline constexpr BlockView& operator=(BlockView const& other) = default;

    /** Move-assigns a new value to this BlockView
     * @param other the BlockView to move
     * @return *this
     */
    inline constexpr BlockView& operator=(BlockView&& other) = default;

    /** Slice out some dimensions
     * @param slices the coordinates to
     */
    template <class... QueryMeshes>
    inline constexpr auto operator[](MCoord<QueryMeshes...> const& slice_spec) const
    {
        auto subview = std::experimental::submdspan(m_raw, get_slicer_for<Meshes>(slice_spec)...);
        using detail::TypeSeq;
        using selected_meshes = type_seq_remove_t<TypeSeq<Meshes...>, TypeSeq<QueryMeshes...>>;
        return ::BlockView(select_by_type_seq<selected_meshes>(m_domain), subview);
    }

    /** Slice out some dimensions
     * @param slices the coordinates to
     */
    template <class... QueryMeshes>
    inline constexpr auto operator[](ProductMDomain<QueryMeshes...> const& odomain) const
    {
        auto subview = std::experimental::submdspan(m_raw, get_slicer_for<Meshes>(odomain)...);
        return ::BlockView(m_domain.restrict(odomain), subview);
    }

    template <class... OMeshes>
    inline constexpr reference operator()(
            TaggedVector<std::size_t, OMeshes> const&... mcoords) const noexcept
    {
        assert(((mcoords >= front<OMeshes>(m_domain)) && ...));
        return m_raw(take_first<Meshes>(mcoords...)...);
    }

    inline constexpr reference operator()(mcoord_type const& indices) const noexcept
    {
        assert(((get<Meshes>(indices) >= front<Meshes>(m_domain)) && ...));
        return m_raw(indices.array());
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
        return m_raw.accessor();
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
                (m_raw.extent(type_seq_rank_v<Meshes, detail::TypeSeq<Meshes...>>)
                 - front<Meshes>(m_domain))...);
    }

    template <class QueryMesh>
    inline constexpr size_type extent() const noexcept
    {
        return m_raw.extent(type_seq_rank_v<QueryMesh, detail::TypeSeq<Meshes...>>)
               - front<QueryMesh>(m_domain);
    }

    inline constexpr size_type size() const noexcept
    {
        return raw_view_without_offset().size();
    }

    inline constexpr size_type unique_size() const noexcept
    {
        return raw_view_without_offset().size();
    }

    static inline constexpr bool is_always_unique() noexcept
    {
        return underlying_raw_view_type::mapping_type::is_always_unique();
    }

    static inline constexpr bool is_always_contiguous() noexcept
    {
        return underlying_raw_view_type::mapping_type::is_always_contiguous();
    }

    static inline constexpr bool is_always_strided() noexcept
    {
        return underlying_raw_view_type::mapping_type::is_always_strided();
    }

    inline constexpr mapping_type mapping() const noexcept
    {
        return m_raw.mapping();
    }

    inline constexpr bool is_unique() const noexcept
    {
        return raw_view_without_offset().is_unique();
    }

    inline constexpr bool is_contiguous() const noexcept
    {
        return raw_view_without_offset().is_contiguous();
    }

    inline constexpr bool is_strided() const noexcept
    {
        return raw_view_without_offset().is_strided();
    }

    template <class QueryMesh>
    inline constexpr auto stride() const
    {
        return m_raw.stride(type_seq_rank_v<QueryMesh, detail::TypeSeq<Meshes...>>);
    }

    /** Swaps this field with another
     * @param other the Block to swap with this one
     */
    inline constexpr void swap(BlockView& other)
    {
        BlockView tmp = std::move(other);
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
        namespace stdex = std::experimental;
        stdex::dextents<mesh_type::rank()> extents_r(::extents<Meshes>(m_domain)...);
        stdex::layout_right::mapping mapping_r(extents_r);

        stdex::dextents<mesh_type::rank()> extents_s(
                (front<Meshes>(m_domain) + ::extents<Meshes>(m_domain))...);
        std::array<std::size_t, mesh_type::rank()> strides_s {
                mapping_r.stride(type_seq_rank_v<Meshes, detail::TypeSeq<Meshes...>>)...};
        stdex::layout_stride::mapping mapping_s(extents_s, strides_s);
        return m_raw.data() + mapping_s(front<Meshes>(m_domain)...);
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr raw_view_type raw_view() const
    {
        return m_raw;
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr underlying_raw_view_type raw_view_without_offset() const
    {
        return underlying_raw_view_type(data(), ::extents<Meshes>(m_domain)...);
    }
};

template <class... Meshes, class ElementType, class Extents, class Layout>
BlockView(
        ProductMDomain<Meshes...> domain,
        std::experimental::mdspan<ElementType, Extents, Layout> raw_view)
        -> BlockView<ProductMDomain<Meshes...>, ElementType, std::experimental::layout_stride>;
