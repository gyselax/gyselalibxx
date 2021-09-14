#pragma once

#include "ddc/block_span.h"

template <class, class>
class Block;

template <class... Meshes, class ElementType>
class Block<ProductMDomain<Meshes...>, ElementType>
    : public BlockSpan<ProductMDomain<Meshes...>, ElementType>
{
public:
    /// type of a span of this full block
    using block_span_type = BlockSpan<ProductMDomain<Meshes...>, ElementType>;

protected:
    /// ND memory view
    using raw_mdspan_type = typename block_span_type::raw_mdspan_type;

public:
    /// type of a view of this full block
    using block_view_type = BlockSpan<ProductMDomain<Meshes...>, ElementType const>;

    /// The dereferenceable part of the co-domain but with indexing starting at 0
    using allocation_mdspan_type = typename block_span_type::allocation_mdspan_type;

    using mdomain_type = typename block_span_type::mdomain_type;

    using mesh_type = ProductMesh<Meshes...>;

    using mcoord_type = typename mdomain_type::mcoord_type;

    using extents_type = typename block_span_type::extents_type;

    using layout_type = typename block_span_type::layout_type;

    using mapping_type = typename block_span_type::mapping_type;

    using element_type = typename block_span_type::element_type;

    using value_type = typename block_span_type::value_type;

    using size_type = typename block_span_type::size_type;

    using difference_type = typename block_span_type::difference_type;

    using pointer = typename block_span_type::pointer;

    using reference = typename block_span_type::reference;

public:
    /** Construct a Block on a domain with uninitialized values
     */
    explicit inline constexpr Block(mdomain_type const& domain)
        : block_span_type(domain, new (std::align_val_t(64)) value_type[domain.size()])
    {
    }

    /** Constructs a new Block by copy
     * 
     * This is deleted, one should use deepcopy
     * @param other the Block to copy
     */
    inline constexpr Block(Block const& other) = delete;

    /** Constructs a new Block by move
     * @param other the Block to move
     */
    inline constexpr Block(Block&& other) = default;

    inline ~Block()
    {
        if (this->m_raw.data()) {
            operator delete(this->data(), std::align_val_t(64));
        }
    }

    /** Copy-assigns a new value to this field
     * @param other the Block to copy
     * @return *this
     */
    inline constexpr Block& operator=(Block const& other) = default;

    /** Move-assigns a new value to this field
     * @param other the Block to move
     * @return *this
     */
    inline constexpr Block& operator=(Block&& other) = default;

    /** Copy-assigns a new value to this field
     * @param other the Block to copy
     * @return *this
     */
    template <class... OMeshes, class OElementType>
    inline Block& operator=(Block<ProductMDomain<OMeshes...>, OElementType>&& other)
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

    template <class... OMeshes>
    inline constexpr element_type const& operator()(
            TaggedVector<std::size_t, OMeshes> const&... mcoords) const noexcept
    {
        return this->m_raw(take_first<Meshes>(mcoords...)...);
    }

    template <class... OMeshes>
    inline constexpr element_type& operator()(
            TaggedVector<std::size_t, OMeshes> const&... mcoords) noexcept
    {
        assert(((mcoords >= front<OMeshes>(this->m_domain)) && ...));
        return this->m_raw(take_first<Meshes>(mcoords...)...);
    }

    inline constexpr element_type const& operator()(mcoord_type const& indices) const noexcept
    {
        assert(((get<Meshes>(indices) >= front<Meshes>(this->m_domain)) && ...));
        return this->m_raw(indices.array());
    }

    inline constexpr element_type& operator()(mcoord_type const& indices) noexcept
    {
        assert(((get<Meshes>(indices) >= front<Meshes>(this->m_domain)) && ...));
        return this->m_raw(indices.array());
    }

    inline constexpr block_view_type cview() const
    {
        return *this;
    }

    inline constexpr block_view_type view() const
    {
        return *this;
    }

    inline constexpr block_span_type view()
    {
        return *this;
    }
};
