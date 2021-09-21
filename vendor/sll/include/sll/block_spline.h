#pragma once

#include <ddc/BlockSpan>
#include <ddc/MDomain>

#include "sll/bsplines.h"

template <class ElementType, class Mesh, std::size_t D>
class Block<ElementType, BSplines<Mesh, D>> : public BlockSpan<ElementType, BSplines<Mesh, D>>
{
public:
    /// ND view on this block
    using block_view_type = BlockSpan<ElementType, BSplines<Mesh, D>>;

    using block_span_type = BlockSpan<ElementType const, BSplines<Mesh, D>>;

    /// ND memory view
    using raw_view_type = std::experimental::mdspan<ElementType, std::experimental::dextents<1>>;

    using bsplines_type = BSplines<Mesh, D>;

    using mcoord_type = typename bsplines_type::mcoord_type;

    using extents_type = typename block_view_type::extents_type;

    using layout_type = typename block_view_type::layout_type;

    using mapping_type = typename block_view_type::mapping_type;

    using element_type = typename block_view_type::element_type;

    using value_type = typename block_view_type::value_type;

    using size_type = typename block_view_type::size_type;

    using difference_type = typename block_view_type::difference_type;

    using pointer = typename block_view_type::pointer;

    using reference = typename block_view_type::reference;

public:
    /** Construct a Block on a domain with uninitialized values
     */
    explicit inline constexpr Block(BSplines<Mesh, D> const& bsplines)
        : block_view_type(
                bsplines,
                raw_view_type(
                        new (std::align_val_t(64)) value_type[bsplines.size()],
                        bsplines.size()))
    {
    }

    /** Constructs a new Block by copy
     * 
     * This is deleted, one should use deepcopy
     * @param other the Block to copy
     */
    inline constexpr Block(const Block& other) = delete;

    /** Constructs a new Block by move
     * @param other the Block to move
     */
    inline constexpr Block(Block&& other) = default;

    inline ~Block()
    {
        if (this->raw_view().data()) {
            operator delete(this->raw_view().data(), std::align_val_t(64));
        }
    }

    /** Copy-assigns a new value to this field
     * @param other the Block to copy
     * @return *this
     */
    inline constexpr Block& operator=(const Block& other) = default;

    /** Move-assigns a new value to this field
     * @param other the Block to move
     * @return *this
     */
    inline constexpr Block& operator=(Block&& other) = default;

    /** Swaps this field with another
     * @param other the Block to swap with this one
     */
    inline constexpr void swap(Block& other)
    {
        Block tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    bsplines_type const& bsplines() const noexcept
    {
        return this->m_bsplines;
    }

    inline constexpr ElementType const& operator()(const mcoord_type& indices) const noexcept
    {
        return this->m_raw(indices);
    }

    inline constexpr ElementType& operator()(const mcoord_type& indices) noexcept
    {
        return this->m_raw(indices);
    }

    inline constexpr block_span_type cview() const
    {
        return *this;
    }

    inline constexpr block_span_type cview()
    {
        return *this;
    }

    inline constexpr block_span_type view() const
    {
        return *this;
    }

    inline constexpr block_view_type view()
    {
        return *this;
    }
};
