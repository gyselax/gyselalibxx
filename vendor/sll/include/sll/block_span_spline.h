#pragma once

#include <ddc/block_span.h>

#include "sll/bsplines.h"

template <class Mesh, std::size_t D, class ElementType>
class BlockSpan<BSplines<Mesh, D>, ElementType>
{
public:
    /// ND memory view
    using raw_view_type = std::experimental::mdspan<ElementType, std::experimental::dextents<1>>;

    using bsplines_type = BSplines<Mesh, D>;

    using mcoord_type = typename bsplines_type::mcoord_type;

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

protected:
    /// The raw view of the data
    raw_view_type m_raw;

    /// The mesh on which this block is defined
    bsplines_type const& m_bsplines;

public:
    /** Constructs a new BlockSpan by copy, yields a new view to the same data
     * @param other the BlockSpan to copy
     */
    inline constexpr BlockSpan(const BlockSpan& other) = default;

    /** Constructs a new BlockSpan by move
     * @param other the BlockSpan to move
     */
    inline constexpr BlockSpan(BlockSpan&& other) = default;

    /** Constructs a new BlockSpan by copy of a block, yields a new view to the same data
     * @param other the BlockSpan to move
     */
    template <class OElementType>
    inline constexpr BlockSpan(BlockSpan<bsplines_type, OElementType> const& other) noexcept
        : m_raw(other.raw_view())
        , m_bsplines(other.bsplines())
    {
    }

    /** Constructs a new BlockSpan from scratch
     * @param bsplines the bsplines that sustains the view
     * @param raw_view the raw view to the data
     */
    inline constexpr BlockSpan(bsplines_type const& bsplines, raw_view_type raw_view)
        : m_raw(raw_view)
        , m_bsplines(bsplines)
    {
    }

    /** Copy-assigns a new value to this BlockSpan, yields a new view to the same data
     * @param other the BlockSpan to copy
     * @return *this
     */
    inline constexpr BlockSpan& operator=(const BlockSpan& other) = default;

    /** Move-assigns a new value to this BlockSpan
     * @param other the BlockSpan to move
     * @return *this
     */
    inline constexpr BlockSpan& operator=(BlockSpan&& other) = default;

    inline constexpr reference operator()(const mcoord_type& indices) const noexcept
    {
        return m_raw(indices);
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
    inline constexpr void swap(BlockSpan& other)
    {
        BlockSpan tmp = std::move(other);
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
