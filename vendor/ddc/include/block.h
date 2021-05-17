#pragma once

#include "blockview.h"

template <class... Tags, class ElementType>
class Block<MDomain<Tags...>, ElementType> : public BlockView<MDomain<Tags...>, ElementType>
{
public:
    /// ND view on this block
    using BlockView_ = BlockView<MDomain<Tags...>, ElementType>;

    /// ND memory view
    using RawView = SpanND<sizeof...(Tags), ElementType>;

    using MDomain_ = MDomain<Tags...>;

    using Mesh = typename MDomain_::Mesh_;

    using MCoord_ = typename MDomain_::MCoord_;

    using extents_type = typename BlockView_::extents_type;

    using layout_type = typename BlockView_::layout_type;

    using accessor_type = typename BlockView_::accessor_type;

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
    template <class... OTags>
    explicit inline constexpr Block(const MDomain<OTags...>& domain)
        : BlockView_(
                domain.mesh(),
                RawView(new (std::align_val_t(64)) value_type[domain.size()],
                        ExtentsND<sizeof...(Tags)>(domain.template extent<Tags>()...)))
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
        delete[] this->raw_view().data();
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
    template <class... OTags, class OElementType>
    inline Block& operator=(Block<MDomain<OTags...>, OElementType>&& other)
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
};

template <class ElementType>
using BlockX = Block<MDomain<Dim::X>, ElementType>;

using DBlockX = BlockX<double>;

template <class ElementType>
using BlockVx = Block<MDomain<Dim::Vx>, ElementType>;

using DBlockVx = BlockVx<double>;

template <class ElementType>
using BlockXVx = Block<MDomain<Dim::X, Dim::Vx>, ElementType>;

using DBlockXVx = BlockXVx<double>;
