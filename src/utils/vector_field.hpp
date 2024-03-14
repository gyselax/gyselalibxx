// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "vector_field_common.hpp"

/**
 * @brief Pre-declaration of VectorField.
 */
template <class ElementType, class Domain, class, class Allocator = ddc::HostAllocator<ElementType>>
class VectorField;


template <class ElementType, class Domain, class DimSeq, class Allocator>
inline constexpr bool enable_field<VectorField<ElementType, Domain, DimSeq, Allocator>> = true;

/**
 * @brief Pre-declaration of VectorFieldSpan.
 */
template <class, class, class, class, class>
class VectorFieldSpan;

/**
 * @brief A class which describes a vector field.
 *
 * A class which describes a vector field. In other words a class which maps a position on a domain
 * to a vector (x,y,z,...). This is done by storing the values at the positions in individual
 * Chunks.
 */
template <class ElementType, class Domain, class NDTag, class Allocator>
class VectorField : public VectorFieldCommon<ddc::Chunk<ElementType, Domain, Allocator>, NDTag>
{
public:
    /**
     * @brief Type describing the object which can be extracted from this VectorField using the get<> function.
     */
    using chunk_type = ddc::Chunk<ElementType, Domain, Allocator>;

private:
    using base_type = VectorFieldCommon<chunk_type, NDTag>;

public:
    /// The type of an element in one of the Chunks comprising the VectorField
    using element_type = typename base_type::element_type;

public:
    /**
     * @brief A type which can hold a reference to this VectorField.
     */
    using span_type = VectorFieldSpan<
            ElementType,
            Domain,
            NDTag,
            std::experimental::layout_right,
            typename Allocator::memory_space>;

    /**
     * @brief A type which can hold a constant reference to this VectorField.
     */
    using view_type = VectorFieldSpan<
            const ElementType,
            Domain,
            NDTag,
            std::experimental::layout_right,
            typename Allocator::memory_space>;

    /**
     * @brief The type of the domain on which the field is defined.
     */
    using mdomain_type = typename base_type::mdomain_type;

    /**
     * @brief The type of the memory space where the field is saved (CPU vs GPU).
     */
    using memory_space = typename chunk_type::memory_space;

private:
    /// Construct a VectorField on a domain with uninitialized values
    template <std::size_t... Is>
    explicit VectorField(
            mdomain_type const& domain,
            Allocator allocator,
            std::index_sequence<Is...> const&)
        : base_type(((void)Is, chunk_type(domain, allocator))...)
    {
    }

    /**
     * Construct a VectorField from a deepcopy of a VectorFieldSpan
     *
     * @param[inout] field_span
     */
    template <class OElementType, class LayoutDomain, class MemorySpace, std::size_t... Is>
    explicit VectorField(
            VectorFieldSpan<OElementType, Domain, NDTag, LayoutDomain, MemorySpace> field_span,
            std::index_sequence<Is...> const&)
        : base_type(chunk_type(ddcHelper::get<ddc::type_seq_element_t<Is, NDTag>>(field_span))...)
    {
    }

    /** Element access using a multi-dimensional DiscreteElement
     * @param delems discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims, typename T, T... ints>
    element_type operator()(
            ddc::DiscreteElement<ODDims...> const& delems,
            std::integer_sequence<T, ints...>) const noexcept
    {
        return element_type((base_type::m_values[ints](delems))...);
    }

public:
    /// Empty VectorField
    VectorField() = default;

    /**
     * Construct a VectorField on a domain with uninitialized values
     *
     * @param[in] domain The domain on which the chunk will be defined.
     * @param[in] allocator An optional allocator used to create the chunks.
     */
    explicit VectorField(mdomain_type const& domain, Allocator allocator = Allocator())
        : VectorField(domain, allocator, std::make_index_sequence<base_type::NDims> {})
    {
    }

    /**
     * Construct a VectorField from a deepcopy of a VectorFieldSpan
     *
     * @param[inout] field_span A reference to another vector field.
     */
    template <class OElementType, class LayoutDomain, class MemorySpace>
    explicit VectorField(
            VectorFieldSpan<OElementType, Domain, NDTag, LayoutDomain, MemorySpace> field_span)
        : VectorField(field_span, std::make_index_sequence<base_type::NDims> {})
    {
    }

    /// Deleted: use deepcopy instead
    VectorField(VectorField const& other) = delete;

    /**
     * Constructs a new VectorField by move
     * @param other the VectorField to move
     */
    VectorField(VectorField&& other) = default;

    /**
     * Copy-assigns a new value to this VectorFieldSpan, yields a new view to the same data
     * @param other the VectorFieldSpan to copy
     * @return *this
     */
    VectorField& operator=(VectorField const& other) = default;

    /**
     * Move-assigns a new value to this VectorFieldSpan
     * @param other the VectorFieldSpan to move
     * @return *this
     */
    VectorField& operator=(VectorField&& other) = default;

    ~VectorField() = default;

    /**
     * Get a constant reference to this vector field.
     *
     * @return A constant reference to this vector field.
     */
    view_type span_cview() const
    {
        return view_type(*this);
    }

    /**
     * Get a constant reference to this vector field.
     *
     * @return A constant reference to this vector field.
     */
    view_type span_view() const
    {
        return view_type(*this);
    }

    /**
     * Get a modifiable reference to this vector field.
     *
     * @return A modifiable reference to this vector field.
     */
    span_type span_view()
    {
        return span_type(*this);
    }

    /** Element access using a list of DiscreteElement
     * @param delems 1D discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims>
    element_type operator()(ddc::DiscreteElement<ODDims> const&... delems) const noexcept
    {
        ddc::DiscreteElement<ODDims...> delem_idx(delems...);
        return this->
        operator()(delem_idx, std::make_integer_sequence<int, element_type::size()> {});
    }

    /** Element access using a multi-dimensional DiscreteElement
     * @param delems discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims, class = std::enable_if_t<sizeof...(ODDims) != 1>>
    element_type operator()(ddc::DiscreteElement<ODDims...> const& delems) const noexcept
    {
        return this->operator()(delems, std::make_integer_sequence<int, element_type::size()> {});
    }


    /**
     * @brief Slice out some dimensions.
     *
     * Get the VectorField on the reduced domain which is obtained by indexing
     * the dimensions QueryDDims at the position slice_spec.
     *
     * @param[in] slice_spec The slice describing the domain of interest.
     *
     * @return A constant reference to the vector field on the sliced domain.
     */
    template <class... QueryDDims>
    auto operator[](ddc::DiscreteElement<QueryDDims...> const& slice_spec) const
    {
        return span_cview()[slice_spec];
    }

    /**
     * @brief Slice out some dimensions.
     *
     * Get the VectorField on the reduced domain which is obtained by indexing
     * the dimensions QueryDDims at the position slice_spec.
     *
     * @param[in] slice_spec The slice describing the domain of interest.
     *
     * @return A modifiable reference to the vector field on the sliced domain.
     */
    template <class... QueryDDims>
    auto operator[](ddc::DiscreteElement<QueryDDims...> const& slice_spec)
    {
        return span_view()[slice_spec];
    }

    /**
     * @brief Slice out some dimensions.
     *
     * Get the VectorField on the reduced domain passed as an argument.
     *
     * @param[in] odomain The domain of interest.
     *
     * @return A modifiable reference to the vector field on the sliced domain.
     */
    template <class... QueryDDims>
    auto operator[](ddc::DiscreteDomain<QueryDDims...> const& odomain) const
    {
        return span_cview()[odomain];
    }

    /**
     * @brief Slice out some dimensions.
     *
     * Get the VectorField on the reduced domain passed as an argument.
     *
     * @param[in] odomain The domain of interest.
     *
     * @return A modifiable reference to the vector field on the sliced domain.
     */
    template <class... QueryDDims>
    auto operator[](ddc::DiscreteDomain<QueryDDims...> const& odomain)
    {
        return span_view()[odomain];
    }
};
