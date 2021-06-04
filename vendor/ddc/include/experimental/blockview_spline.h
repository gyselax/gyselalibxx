#pragma once

#include "blockview.h"
#include "bsplines_uniform.h"

template <class Tag, std::size_t D, class ElementType, bool CONTIGUOUS>
class BlockView<experimental::UniformBSplines<Tag, D>, ElementType, CONTIGUOUS>
{
public:
    /// ND memory view
    using raw_view_type = SpanND<1, ElementType, CONTIGUOUS>;

    using bsplines_type = experimental::UniformBSplines<Tag, D>;

    using mcoord_type = typename bsplines_type::mcoord_type;

    using extents_type = typename raw_view_type::extents_type;

    using layout_type = typename raw_view_type::layout_type;

    using accessor_type = typename raw_view_type::accessor_type;

    using mapping_type = typename raw_view_type::mapping_type;

    using element_type = typename raw_view_type::element_type;

    using value_type = typename raw_view_type::value_type;

    using index_type = typename raw_view_type::index_type;

    using difference_type = typename raw_view_type::difference_type;

    using pointer = typename raw_view_type::pointer;

    using reference = typename raw_view_type::reference;

protected:
    /// The raw view of the data
    raw_view_type m_raw;

    /// The mesh on which this block is defined
    bsplines_type m_bsplines;

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
    inline constexpr BlockView(
            BlockView<bsplines_type, OElementType, CONTIGUOUS> const& other) noexcept
        : m_raw(other.raw_view())
        , m_bsplines(other.bsplines())
    {
    }

    /** Constructs a new BlockView from scratch
     * @param bsplines the bsplines that sustains the view
     * @param raw_view the raw view to the data
     */
    inline constexpr BlockView(bsplines_type const& bsplines, raw_view_type raw_view)
        : m_raw(raw_view)
        , m_bsplines(bsplines)
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

    static constexpr std::size_t tag_rank()
    {
        return 0;
    }

    template <class... IndexType>
    inline constexpr reference operator()(IndexType&&... indices) const noexcept
    {
        return m_raw(std::forward<IndexType>(indices)...);
    }

    inline constexpr reference operator()(const mcoord_type& indices) const noexcept
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

    /** Provide access to the bsplines on which this block is defined
     * @return the bsplines on which this block is defined
     */
    inline constexpr bsplines_type bsplines() const noexcept
    {
        return m_bsplines;
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr raw_view_type raw_view()
    {
        return m_raw;
    }

    /** Provide a constant view of the data
     * @return a constant view of the data
     */
    inline constexpr const raw_view_type raw_view() const
    {
        return m_raw;
    }
};
