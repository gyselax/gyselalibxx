#pragma once

#include "ddc/chunk_span.hpp"

template <class, class>
class Chunk;

template <class ElementType, class... DDims>
class Chunk<ElementType, DiscreteDomain<DDims...>>
    : public ChunkSpan<ElementType, DiscreteDomain<DDims...>>
{
public:
    /// type of a span of this full chunk
    using span_type = ChunkSpan<ElementType, DiscreteDomain<DDims...>>;

protected:
    /// ND memory view
    using internal_mdspan_type = typename span_type::internal_mdspan_type;

public:
    /// type of a view of this full chunk
    using view_type = ChunkSpan<ElementType const, DiscreteDomain<DDims...>>;

    /// The dereferenceable part of the co-domain but with indexing starting at 0
    using allocation_mdspan_type = typename span_type::allocation_mdspan_type;

    using mdomain_type = typename span_type::mdomain_type;

    using mcoord_type = typename mdomain_type::mcoord_type;

    using extents_type = typename span_type::extents_type;

    using layout_type = typename span_type::layout_type;

    using mapping_type = typename span_type::mapping_type;

    using element_type = typename span_type::element_type;

    using value_type = typename span_type::value_type;

    using size_type = typename span_type::size_type;

    using difference_type = typename span_type::difference_type;

    using pointer = typename span_type::pointer;

    using reference = typename span_type::reference;

public:
    /** Construct a Chunk on a domain with uninitialized values
     */
    explicit inline constexpr Chunk(mdomain_type const& domain)
        : span_type(new (std::align_val_t(64)) value_type[domain.size()], domain)
    {
    }

    /// Deleted: use deepcopy instead
    inline constexpr Chunk(Chunk const& other) = delete;

    /** Constructs a new Chunk by move
     * @param other the Chunk to move
     */
    inline constexpr Chunk(Chunk&& other) : span_type(std::move(other))
    {
        other.m_internal_mdspan = internal_mdspan_type();
    }

    inline ~Chunk()
    {
        if (this->m_internal_mdspan.data()) {
            operator delete(this->data(), std::align_val_t(64));
        }
    }

    /// Deleted: use deepcopy instead
    inline constexpr Chunk& operator=(Chunk const& other) = delete;

    /** Move-assigns a new value to this field
     * @param other the Chunk to move
     * @return *this
     */
    inline constexpr Chunk& operator=(Chunk&& other)
    {
        static_cast<span_type&>(*this) = std::move(other);
        other.m_internal_mdspan = internal_mdspan_type();
        return *this;
    }

    /** Swaps this field with another
     * @param other the Chunk to swap with this one
     */
    inline constexpr void swap(Chunk& other)
    {
        Chunk tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    // Warning: Do not use DiscreteCoordinate because of template deduction issue with clang 12
    template <class... ODDims>
    inline constexpr element_type const& operator()(
            detail::TaggedVector<DiscreteCoordElement, ODDims> const&... mcoords) const noexcept
    {
        return this->m_internal_mdspan(take<DDims>(mcoords...)...);
    }

    // Warning: Do not use DiscreteCoordinate because of template deduction issue with clang 12
    template <class... ODDims>
    inline constexpr element_type& operator()(
            detail::TaggedVector<DiscreteCoordElement, ODDims> const&... mcoords) noexcept
    {
        assert(((mcoords >= front<ODDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(take<DDims>(mcoords...)...);
    }

    inline constexpr element_type const& operator()(mcoord_type const& indices) const noexcept
    {
        assert(((get<DDims>(indices) >= front<DDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(indices.array());
    }

    inline constexpr element_type& operator()(mcoord_type const& indices) noexcept
    {
        assert(((get<DDims>(indices) >= front<DDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(indices.array());
    }

    inline constexpr view_type span_cview() const
    {
        return *this;
    }

    inline constexpr view_type span_view() const
    {
        return *this;
    }

    inline constexpr span_type span_view()
    {
        return *this;
    }
};
