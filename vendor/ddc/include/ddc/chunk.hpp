// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include "ddc/chunk_common.hpp"
#include "ddc/chunk_span.hpp"

template <class ElementType, class, class Allocator = std::allocator<ElementType>>
class Chunk;

template <class ElementType, class SupportType, class Allocator>
inline constexpr bool enable_chunk<Chunk<ElementType, SupportType, Allocator>> = true;

template <class ElementType, class... DDims, class Allocator>
class Chunk<ElementType, DiscreteDomain<DDims...>, Allocator>
    : public ChunkCommon<ElementType, DiscreteDomain<DDims...>, std::experimental::layout_right>
{
protected:
    using base_type
            = ChunkCommon<ElementType, DiscreteDomain<DDims...>, std::experimental::layout_right>;

    /// ND memory view
    using internal_mdspan_type = typename base_type::internal_mdspan_type;

public:
    /// type of a span of this full chunk
    using span_type
            = ChunkSpan<ElementType, DiscreteDomain<DDims...>, std::experimental::layout_right>;

    /// type of a view of this full chunk
    using view_type = ChunkSpan<
            ElementType const,
            DiscreteDomain<DDims...>,
            std::experimental::layout_right>;

    /// The dereferenceable part of the co-domain but with indexing starting at 0
    using allocation_mdspan_type = typename base_type::allocation_mdspan_type;

    using const_allocation_mdspan_type = typename base_type::const_allocation_mdspan_type;

    using mdomain_type = typename base_type::mdomain_type;

    using mcoord_type = typename base_type::mcoord_type;

    using extents_type = typename base_type::extents_type;

    using layout_type = typename base_type::layout_type;

    using mapping_type = typename base_type::mapping_type;

    using element_type = typename base_type::element_type;

    using value_type = typename base_type::value_type;

    using size_type = typename base_type::size_type;

    using difference_type = typename base_type::difference_type;

    using pointer = typename base_type::pointer;

    using reference = typename base_type::reference;

    template <class, class, class>
    friend class Chunk;

private:
    Allocator m_allocator;

public:
    /// Empty Chunk
    Chunk() = default;

    /// Construct a Chunk on a domain with uninitialized values
    explicit Chunk(mdomain_type const& domain, Allocator allocator = Allocator())
        : base_type(std::allocator_traits<Allocator>::allocate(m_allocator, domain.size()), domain)
        , m_allocator(std::move(allocator))
    {
    }

    /// Construct a Chunk from a deepcopy of a ChunkSpan
    template <class OElementType, class... ODDims, class LayoutType>
    explicit Chunk(
            ChunkSpan<OElementType, DiscreteDomain<ODDims...>, LayoutType> chunk_span,
            Allocator allocator = Allocator())
        : Chunk(chunk_span.domain(), std::move(allocator))
    {
        deepcopy(span_view(), chunk_span);
    }

    /// Deleted: use deepcopy instead
    Chunk(Chunk const& other) = delete;

    /** Constructs a new Chunk by move
     * @param other the Chunk to move
     */
    Chunk(Chunk&& other)
        : base_type(std::move(static_cast<base_type&>(other)))
        , m_allocator(std::move(other.m_allocator))
    {
        other.m_internal_mdspan = internal_mdspan_type();
    }

    ~Chunk()
    {
        if (this->m_internal_mdspan.data()) {
            std::allocator_traits<Allocator>::deallocate(m_allocator, this->data(), this->size());
        }
    }

    /// Deleted: use deepcopy instead
    Chunk& operator=(Chunk const& other) = delete;

    /** Move-assigns a new value to this field
     * @param other the Chunk to move
     * @return *this
     */
    Chunk& operator=(Chunk&& other)
    {
        assert(this != &other);
        if (this->m_internal_mdspan.data()) {
            std::allocator_traits<Allocator>::deallocate(m_allocator, this->data(), this->size());
        }
        static_cast<base_type&>(*this) = std::move(static_cast<base_type&>(other));
        m_allocator = std::move(other.m_allocator);

        other.m_internal_mdspan = internal_mdspan_type();
        return *this;
    }

    /// Slice out some dimensions
    template <class... QueryDDims>
    auto operator[](DiscreteCoordinate<QueryDDims...> const& slice_spec) const
    {
        return view_type(*this)[slice_spec];
    }

    /// Slice out some dimensions
    template <class... QueryDDims>
    auto operator[](DiscreteCoordinate<QueryDDims...> const& slice_spec)
    {
        return span_view()[slice_spec];
    }

    /// Slice out some dimensions
    template <class... QueryDDims>
    auto operator[](DiscreteDomain<QueryDDims...> const& odomain) const
    {
        return span_view()[odomain];
    }

    /// Slice out some dimensions
    template <class... QueryDDims>
    auto operator[](DiscreteDomain<QueryDDims...> const& odomain)
    {
        return span_view()[odomain];
    }

    /** Element access using a list of DiscreteCoordinate
     * @param mcoords 1D discrete coordinates
     * @return const-reference to this element
     */
    template <class... ODDims>
    element_type const& operator()(DiscreteCoordinate<ODDims> const&... mcoords) const noexcept
    {
        static_assert(sizeof...(ODDims) == sizeof...(DDims), "Invalid number of dimensions");
        assert(((mcoords >= front<ODDims>(this->m_domain)) && ...));
        assert(((mcoords <= back<ODDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(take<DDims>(mcoords...).uid()...);
    }

    /** Element access using a list of DiscreteCoordinate
     * @param mcoords 1D discrete coordinates
     * @return reference to this element
     */
    template <class... ODDims>
    element_type& operator()(DiscreteCoordinate<ODDims> const&... mcoords) noexcept
    {
        static_assert(sizeof...(ODDims) == sizeof...(DDims), "Invalid number of dimensions");
        assert(((mcoords >= front<ODDims>(this->m_domain)) && ...));
        assert(((mcoords <= back<ODDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(take<DDims>(mcoords...).uid()...);
    }

    /** Element access using a multi-dimensional DiscreteCoordinate
     * @param mcoord discrete coordinates
     * @return const-reference to this element
     */
    template <class... ODDims, class = std::enable_if_t<sizeof...(ODDims) != 1>>
    element_type const& operator()(DiscreteCoordinate<ODDims...> const& mcoord) const noexcept
    {
        static_assert(sizeof...(ODDims) == sizeof...(DDims), "Invalid number of dimensions");
        assert(((select<ODDims>(mcoord) >= front<ODDims>(this->m_domain)) && ...));
        assert(((select<ODDims>(mcoord) <= back<ODDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(uid<DDims>(mcoord)...);
    }

    /** Element access using a multi-dimensional DiscreteCoordinate
     * @param mcoord discrete coordinates
     * @return reference to this element
     */
    template <class... ODDims, class = std::enable_if_t<sizeof...(ODDims) != 1>>
    element_type& operator()(DiscreteCoordinate<ODDims...> const& mcoord) noexcept
    {
        static_assert(sizeof...(ODDims) == sizeof...(DDims), "Invalid number of dimensions");
        assert(((select<ODDims>(mcoord) >= front<ODDims>(this->m_domain)) && ...));
        assert(((select<ODDims>(mcoord) <= back<ODDims>(this->m_domain)) && ...));
        return this->m_internal_mdspan(uid<DDims>(mcoord)...);
    }

    /** Access to the underlying allocation pointer
     * @return read-only allocation pointer
     */
    ElementType const* data() const
    {
        return base_type::data();
    }

    /** Access to the underlying allocation pointer
     * @return allocation pointer
     */
    ElementType* data()
    {
        return base_type::data();
    }

    /** Provide a mdspan on the memory allocation
     * @return read-only allocation mdspan
     */
    const_allocation_mdspan_type allocation_mdspan() const
    {
        return base_type::allocation_mdspan();
    }

    /** Provide a mdspan on the memory allocation
     * @return allocation mdspan
     */
    allocation_mdspan_type allocation_mdspan()
    {
        return base_type::allocation_mdspan();
    }

    view_type span_cview() const
    {
        return view_type(*this);
    }

    view_type span_view() const
    {
        return view_type(*this);
    }

    span_type span_view()
    {
        return span_type(*this);
    }
};
