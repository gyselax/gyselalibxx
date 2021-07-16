#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include "mdomain.h"
#include "product_mdomain.h"
#include "product_mesh.h"
#include "taggedarray.h"
#include "view.h"

template <class, class>
class Block;

template <class, class, bool = true>
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

template <class... Meshes, class ElementType, bool CONTIGUOUS>
class BlockView<ProductMDomain<Meshes...>, ElementType, CONTIGUOUS>
{
public:
    using mdomain_type = ProductMDomain<Meshes...>;

    using mesh_type = ProductMesh<Meshes...>;

    /// ND memory view
    using raw_view_type = SpanND<mesh_type::rank(), ElementType, CONTIGUOUS>;

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

    template <class, class, bool>
    friend class BlockView;

protected:
    /// This adaptor transforms the spec from `(start, count)` to `[begin, end[`
    template <class SliceSpec>
    static SliceSpec slice_spec_adaptor(SliceSpec const& slice_spec)
    {
        if constexpr (std::is_convertible_v<SliceSpec, std::pair<std::size_t, std::size_t>>) {
            return std::pair(slice_spec.first, slice_spec.first + slice_spec.second);
        } else {
            return slice_spec;
        }
    }

    template <class QueryMesh, class... OMeshes>
    static auto get_slicer_for(MCoord<OMeshes...> const& c)
    {
        if constexpr (has_tag_v<QueryMesh, MCoord<OMeshes...>>) {
            return c.template get<QueryMesh>();
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
        : m_raw(other.raw_view())
        , m_domain(other.domain())
    {
    }

    /** Constructs a new BlockView by copy of a block, yields a new view to the same data
     * @param other the BlockView to move
     */
    template <class OElementType>
    inline constexpr BlockView(
            BlockView<mdomain_type, OElementType, CONTIGUOUS> const& other) noexcept
        : m_raw(other.raw_view())
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
    inline constexpr auto operator[](MCoord<QueryMeshes...> mcoord) const
    {
        return this->subblockview(get_slicer_for<Meshes>(mcoord)...);
    }

    template <
            class... IndexType,
            std::enable_if_t<(... && std::is_convertible_v<IndexType, std::size_t>), int> = 0>
    inline constexpr reference operator()(IndexType&&... indices) const noexcept
    {
        return m_raw(std::forward<IndexType>(indices)...);
    }

    inline constexpr reference operator()(mcoord_type const& indices) const noexcept
    {
        return m_raw(indices.array());
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

    static inline constexpr size_type static_extent(size_t r) noexcept
    {
        return extents_type::static_extent(r);
    }

    inline constexpr extents_type extents() const noexcept
    {
        return m_raw.extents();
    }

    inline constexpr size_type extent(size_t dim) const noexcept
    {
        return m_raw.extent(dim);
    }

    inline constexpr size_type size() const noexcept
    {
        return m_raw.size();
    }

    inline constexpr size_type unique_size() const noexcept
    {
        return m_raw.unique_size();
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
        return m_raw.mapping();
    }

    inline constexpr bool is_unique() const noexcept
    {
        return m_raw.is_unique();
    }

    inline constexpr bool is_contiguous() const noexcept
    {
        return m_raw.is_contiguous();
    }

    inline constexpr bool is_strided() const noexcept
    {
        return m_raw.is_strided();
    }

    inline constexpr size_type stride(size_t r) const
    {
        return m_raw.stride(r);
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

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr raw_view_type raw_view() const
    {
        return m_raw;
    }

    /** Slice out some dimensions
     * @param slices the coordinates to 
     */
    template <class... SliceSpecs>
    inline constexpr auto subblockview(SliceSpecs&&... slices) const
    {
        static_assert(sizeof...(SliceSpecs) == sizeof...(Meshes));
        auto subview = std::experimental::submdspan(m_raw, slice_spec_adaptor(slices)...);
        return ::BlockView(m_domain.subdomain(slices...), subview);
    }
};

template <class... Meshes, class ElementType, class Extents, class LayoutPolicy>
BlockView(
        ProductMDomain<Meshes...> domain,
        std::experimental::mdspan<ElementType, Extents, LayoutPolicy> raw_view)
        -> BlockView<
                ProductMDomain<Meshes...>,
                ElementType,
                std::is_same_v<LayoutPolicy, std::experimental::layout_right>>;
