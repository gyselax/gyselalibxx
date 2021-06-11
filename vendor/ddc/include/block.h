#pragma once

#include "blockview.h"

template <class, class>
class Block;

template <class Mesh, class ElementType>
class Block<MDomainImpl<Mesh>, ElementType> : public BlockView<MDomainImpl<Mesh>, ElementType>
{
public:
    /// ND view on this block
    using BlockView_ = BlockView<MDomainImpl<Mesh>, ElementType>;

    using BlockSpan_ = BlockView<MDomainImpl<Mesh>, ElementType const>;

    /// ND memory view
    using RawView = typename BlockView_::RawView;

    using MDomain_ = typename BlockView_::MDomain_;

    using Mesh_ = Mesh;

    using MCoord_ = typename MDomain_::MCoord_;

    using extents_type = typename BlockView_::extents_type;

    using layout_type = typename BlockView_::layout_type;

    using mapping_type = typename BlockView_::mapping_type;

    using element_type = typename BlockView_::element_type;

    using value_type = typename BlockView_::value_type;

    using index_type = typename BlockView_::index_type;

    using difference_type = typename BlockView_::difference_type;

    using pointer = typename BlockView_::pointer;

    using reference = typename BlockView_::reference;

    template <class, class, bool>
    friend class BlockView;

public:
    /** Construct a Block on a domain with uninitialized values
     */
    explicit inline constexpr Block(const MDomain_& domain)
        : BlockView_(
                domain.mesh(),
                RawView(new (std::align_val_t(64)) value_type[domain.size()], domain.extents()))
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

    /** Copy-assigns a new value to this field
     * @param other the Block to copy
     * @return *this
     */
    template <class OMesh, class OElementType>
    inline Block& operator=(Block<MDomainImpl<OMesh>, OElementType>&& other)
    {
        copy(*this, other);
        return *this;
    }

    /** Swaps this field with another
     * @param other the Block to swap with this one
     */
    inline constexpr void swap(Block& other)
    {
        Block tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    template <class... IndexType>
    inline constexpr ElementType const& operator()(IndexType&&... indices) const noexcept
    {
        return this->m_raw(std::forward<IndexType>(indices)...);
    }

    template <class... IndexType>
    inline constexpr ElementType& operator()(IndexType&&... indices) noexcept
    {
        return this->m_raw(std::forward<IndexType>(indices)...);
    }

    inline constexpr ElementType const& operator()(const MCoord_& indices) const noexcept
    {
        return this->m_raw(indices.array());
    }

    inline constexpr ElementType& operator()(const MCoord_& indices) noexcept
    {
        return this->m_raw(indices.array());
    }

    inline constexpr BlockView_ cview() const
    {
        return *this;
    }

    inline constexpr BlockView_ cview()
    {
        return *this;
    }

    inline constexpr BlockView_ view() const
    {
        return *this;
    }

    inline constexpr BlockSpan_ view()
    {
        return *this;
    }
};
