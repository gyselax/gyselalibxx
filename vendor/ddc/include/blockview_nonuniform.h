#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include "blockview.h"
#include "mdomain.h"
#include "nonuniformmesh.h"
#include "taggedarray.h"
#include "view.h"

template <class... Tags, class ElementType, bool CONTIGUOUS>
class BlockView<MDomainImpl<NonUniformMesh<Tags...>>, ElementType, CONTIGUOUS>
{
public:
    /// ND memory view
    using RawView = SpanND<sizeof...(Tags), ElementType, CONTIGUOUS>;

    using MDomain_ = MDomainImpl<NonUniformMesh<Tags...>>;

    using Mesh = typename MDomain_::Mesh_;

    using MCoord_ = typename MDomain_::MCoord_;

    using extents_type = typename RawView::extents_type;

    using layout_type = typename RawView::layout_type;

    using accessor_type = typename RawView::accessor_type;

    using mapping_type = typename RawView::mapping_type;

    using element_type = typename RawView::element_type;

    using value_type = typename RawView::value_type;

    using index_type = typename RawView::index_type;

    using difference_type = typename RawView::difference_type;

    using pointer = typename RawView::pointer;

    using reference = typename RawView::reference;

    template <class, class, bool>
    friend class BlockView;

protected:
    /// The raw view of the data
    RawView m_raw;

    /// The mesh on which this block is defined
    Mesh m_mesh;

public:
    /** Constructs a new BlockView by copy, yields a new view to the same data
     * @param other the BlockView to copy
     */
    inline constexpr BlockView(const BlockView& other) = default;

    /** Constructs a new BlockView by move
     * @param other the BlockView to move
     */
    inline constexpr BlockView(BlockView&& other) = default;

    /** Constructs a new BlockView by copy of a block, yields a new view to the same data
     * @param other the BlockView to move
     */
    template <class OElementType>
    inline constexpr BlockView(const Block<MDomain_, OElementType>& other) noexcept
        : m_raw(other.raw_view())
        , m_mesh(other.mesh())
    {
    }

    /** Constructs a new BlockView by copy of a block, yields a new view to the same data
     * @param other the BlockView to move
     */
    template <class OElementType>
    inline constexpr BlockView(const BlockView<MDomain_, OElementType, CONTIGUOUS>& other) noexcept
        : m_raw(other.raw_view())
        , m_mesh(other.mesh())
    {
    }

    /** Constructs a new BlockView from scratch
     * @param mesh the mesh that sustains the view
     * @param raw_view the raw view to the data
     */
    inline constexpr BlockView(const Mesh& mesh, RawView raw_view) : m_raw(raw_view), m_mesh(mesh)
    {
    }

    /** Copy-assigns a new value to this BlockView, yields a new view to the same data
     * @param other the BlockView to copy
     * @return *this
     */
    inline constexpr BlockView& operator=(const BlockView& other) = default;

    /** Move-assigns a new value to this BlockView
     * @param other the BlockView to move
     * @return *this
     */
    inline constexpr BlockView& operator=(BlockView&& other) = default;

    template <class QueryTag>
    static constexpr std::size_t tag_rank()
    {
        return detail::RankIn<detail::SingleType<QueryTag>, detail::TypeSeq<Tags...>>::val;
    }

    template <class... IndexType>
    inline constexpr reference operator()(IndexType&&... indices) const noexcept
    {
        return m_raw(std::forward<IndexType>(indices)...);
    }

    template <class... OTags>
    inline constexpr reference operator()(const MCoord<OTags...>& indices) const noexcept
    {
        return m_raw(MCoord_(indices).array());
    }

    inline constexpr reference operator()(const MCoord_& indices) const noexcept
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

    static inline constexpr index_type static_extent(size_t r) noexcept
    {
        return extents_type::static_extent(r);
    }

    inline constexpr extents_type extents() const noexcept
    {
        return m_raw.extents();
    }

    inline constexpr index_type extent(size_t dim) const noexcept
    {
        return m_raw.extent(dim);
    }

    inline constexpr index_type size() const noexcept
    {
        return m_raw.size();
    }

    inline constexpr index_type unique_size() const noexcept
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

    inline constexpr index_type stride(size_t r) const
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
    inline constexpr Mesh mesh() const noexcept
    {
        return m_mesh;
    }

    /** Provide access to the domain on which this block is defined
     * @return the domain on which this block is defined
     */
    inline constexpr MDomain_ domain() const noexcept
    {
        return MDomain_(mesh(), mcoord_end<Tags...>(raw_view().extents()));
    }

    /** Provide access to the domain on which this block is defined
     * @return the domain on which this block is defined
     */
    template <class... OTags>
    inline constexpr MDomainImpl<NonUniformMesh<OTags...>> domain() const noexcept
    {
        return MDomainImpl<NonUniformMesh<OTags...>>(domain());
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr RawView raw_view()
    {
        return m_raw;
    }

    /** Provide a constant view of the data
     * @return a constant view of the data
     */
    inline constexpr const RawView raw_view() const
    {
        return m_raw;
    }
};
