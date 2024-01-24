// SPDX-License-Identifier: MIT

#pragma once

#include "vector_field.hpp"


template <
        class ElementType,
        class Domain,
        class NDTag,
        class LayoutStridedPolicy = std::experimental::layout_right,
        class MemorySpace = Kokkos::DefaultHostExecutionSpace::memory_space>
class VectorFieldSpan;

template <
        class ElementType,
        class Domain,
        class NDTag,
        class LayoutStridedPolicy,
        class MemorySpace>
inline constexpr bool enable_field<
        VectorFieldSpan<ElementType, Domain, NDTag, LayoutStridedPolicy, MemorySpace>> = true;

template <
        class ElementType,
        class Domain,
        class NDTag,
        class LayoutStridedPolicy,
        class MemorySpace>
inline constexpr bool enable_borrowed_field<
        VectorFieldSpan<ElementType, Domain, NDTag, LayoutStridedPolicy, MemorySpace>> = true;


template <
        class ElementType,
        class Domain,
        class NDTag,
        class LayoutStridedPolicy,
        class MemorySpace>
class VectorFieldSpan
    : public VectorFieldCommon<
              ddc::ChunkSpan<ElementType, Domain, LayoutStridedPolicy, MemorySpace>,
              NDTag>
{
public:
    /**
     * @brief Type describing the object which can be extracted from this VectorFieldSpan using the get<> function.
     */
    using chunk_type = ddc::ChunkSpan<ElementType, Domain, LayoutStridedPolicy, MemorySpace>;

private:
    using base_type = VectorFieldCommon<chunk_type, NDTag>;

public:
    /// The type of an element in one of the ChunkSpans comprising the VectorFieldSpan
    using element_type = typename base_type::element_type;

public:
    /**
     * @brief A type which can hold a modifiable reference to a VectorField.
     *
     * A type which can hold a reference to a VectorField. If this object is modifiable then
     * so is the span type.
     */
    using span_type = VectorFieldSpan<ElementType, Domain, NDTag, LayoutStridedPolicy, MemorySpace>;
    /**
     * @brief A type which can hold a constant reference to a VectorField.
     */
    using view_type
            = VectorFieldSpan<const ElementType, Domain, NDTag, LayoutStridedPolicy, MemorySpace>;

    /**
     * @brief Type describing the way in which the data is laid out in the Chunk span memory.
     *
     * Type describing the way in which the data is laid out in the Chunk span memory.
     * I.e. it describes whether it is contiguous or not.
     */
    using layout_type = LayoutStridedPolicy;

    /**
     * @brief The type of the domain on which the field is defined.
     */
    using mdomain_type = typename base_type::mdomain_type;

    /**
     * @brief The type of the memory space where the field is saved (CPU vs GPU).
     */
    using memory_space = typename chunk_type::memory_space;

private:
    template <class, class, class, class, class>
    friend class VectorFieldSpan;

    template <class OElementType, class Allocator, std::size_t... Is>
    KOKKOS_FUNCTION constexpr VectorFieldSpan(
            VectorField<OElementType, Domain, NDTag, Allocator>& other,
            std::index_sequence<Is...> const&) noexcept
        : base_type((chunk_type(ddcHelper::get<ddc::type_seq_element_t<Is, NDTag>>(other)))...)
    {
    }

    /** Constructs a new VectorFieldSpan from a VectorField, yields a new view to the same data
     * @param other the VectorField to view
     */
    // Disabled by SFINAE in the case of `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            std::size_t... Is,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            class Allocator>
    KOKKOS_FUNCTION constexpr VectorFieldSpan(
            VectorField<OElementType, Domain, NDTag, Allocator> const& other,
            std::index_sequence<Is...> const&) noexcept
        : base_type((chunk_type(ddcHelper::get<ddc::type_seq_element_t<Is, NDTag>>(other)))...)
    {
    }

    template <
            class... Chunks,
            class = std::enable_if_t<std::conjunction_v<std::is_same<Chunks, chunk_type>...>>>
    KOKKOS_FUNCTION constexpr VectorFieldSpan(Chunks&&... chunks) : base_type(std::move(chunks)...)
    {
    }

    /** Constructs a new VectorFieldSpan by copy of a chunk, yields a new view to the same data
     * @param other the VectorFieldSpan to move
     */
    template <class OElementType, std::size_t... Is>
    KOKKOS_FUNCTION constexpr VectorFieldSpan(
            VectorFieldSpan<
                    OElementType,
                    mdomain_type,
                    NDTag,
                    LayoutStridedPolicy,
                    MemorySpace> const& other,
            std::index_sequence<Is...> const&) noexcept
        : base_type((ddcHelper::get<ddc::type_seq_element_t<Is, NDTag>>(other))...)
    {
    }

    template <class SliceType, std::size_t... Is>
    constexpr auto get_slice(SliceType const& slice_spec, std::index_sequence<Is...> const&)
    {
        auto chunk_slices = std::make_tuple(
                this->template get<ddc::type_seq_element_t<Is, NDTag>>()[slice_spec]...);
        using ChunkType = std::tuple_element_t<0, decltype(chunk_slices)>;
        return VectorFieldSpan<
                ElementType,
                typename ChunkType::mdomain_type,
                NDTag,
                typename ChunkType::layout_type,
                typename ChunkType::memory_space>(std::move(std::get<Is>(chunk_slices))...);
    }

    /** Element access using a multi-dimensional DiscreteElement
     * @param delems discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims, typename T, T... ints>
    KOKKOS_FUNCTION element_type operator()(
            ddc::DiscreteElement<ODDims...> const& delems,
            std::integer_sequence<T, ints...>) const noexcept
    {
        return element_type((base_type::m_values[ints](delems))...);
    }

public:
    /// Empty VectorFieldSpan
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorFieldSpan() = default;

    /// VectorFieldSpan destructor
    KOKKOS_DEFAULTED_FUNCTION ~VectorFieldSpan() = default;

    /** Constructs a new VectorFieldSpan by copy, yields a new view to the same data
     * @param other the VectorFieldSpan to copy
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorFieldSpan(VectorFieldSpan const& other) = default;

    /** Constructs a new VectorFieldSpan by move
     * @param other the VectorFieldSpan to move
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorFieldSpan(VectorFieldSpan&& other) = default;

    /** Constructs a new VectorFieldSpan from a VectorField, yields a new view to the same data
     * @param other the VectorField to view
     */
    template <class OElementType, class Allocator>
    KOKKOS_FUNCTION constexpr VectorFieldSpan(
            VectorField<OElementType, Domain, NDTag, Allocator>& other) noexcept
        : VectorFieldSpan(other, std::make_index_sequence<base_type::NDims> {})
    {
    }

    /** Constructs a new VectorFieldSpan from a VectorField, yields a new view to the same data
     * @param other the VectorField to view
     */
    // Disabled by SFINAE in the case of `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            class Allocator>
    KOKKOS_FUNCTION constexpr VectorFieldSpan(
            VectorField<OElementType, Domain, NDTag, Allocator> const& other) noexcept
        : VectorFieldSpan(other, std::make_index_sequence<base_type::NDims> {})
    {
    }

    /** Constructs a new VectorFieldSpan by copy of a chunk, yields a new view to the same data
     * @param other the VectorFieldSpan to move
     */
    template <class OElementType>
    KOKKOS_FUNCTION constexpr VectorFieldSpan(VectorFieldSpan<
                                              OElementType,
                                              mdomain_type,
                                              NDTag,
                                              LayoutStridedPolicy,
                                              MemorySpace> const& other) noexcept
        : VectorFieldSpan(other, std::make_index_sequence<base_type::NDims> {})
    {
    }

    /** Constructs a new VectorFieldSpan from scratch
     * @param ptr the allocation pointer to the data
     * @param domain the domain that sustains the view
     */
    template <
            class... OElementType,
            class = std::enable_if_t<
                    std::conjunction_v<std::is_same<OElementType, ElementType>...>>,
            class = std::enable_if_t<sizeof...(OElementType) == base_type::NDims>>
    KOKKOS_FUNCTION VectorFieldSpan(mdomain_type const& domain, OElementType*... ptr)
        : base_type((chunk_type(ptr, domain))...)
    {
    }

    /** Copy-assigns a new value to this VectorFieldSpan, yields a new view to the same data
     * @param other the VectorFieldSpan to copy
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorFieldSpan& operator=(VectorFieldSpan const& other)
            = default;

    /** Move-assigns a new value to this VectorFieldSpan
     * @param other the VectorFieldSpan to move
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorFieldSpan& operator=(VectorFieldSpan&& other)
            = default;

    /**
     * Get a constant reference to the vector field referred to by this vector field span.
     *
     * @return A constant reference to the vector field.
     */
    constexpr view_type span_cview() const
    {
        return view_type(*this);
    }

    /**
     * Get a modifiable reference to the vector field referred to by this vector field span.
     *
     * @return A constant reference to the vector field.
     */
    constexpr span_type span_view() const
    {
        return *this;
    }

    /** Element access using a list of DiscreteElement
     * @param delems 1D discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims>
    KOKKOS_FUNCTION element_type
    operator()(ddc::DiscreteElement<ODDims> const&... delems) const noexcept
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
    KOKKOS_FUNCTION element_type
    operator()(ddc::DiscreteElement<ODDims...> const& delems) const noexcept
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
     * @return A reference to the vector field on the sliced domain.
     */
    template <class... QueryDDims>
    constexpr auto operator[](ddc::DiscreteElement<QueryDDims...> const& slice_spec)
    {
        return get_slice(slice_spec, std::make_index_sequence<base_type::NDims> {});
    }

    /**
     * @brief Slice out some dimensions.
     *
     * Get the VectorField on the reduced domain passed as an argument.
     *
     * @param[in] odomain The domain of interest.
     *
     * @return A reference to the vector field on the sliced domain.
     */
    template <class... QueryDDims>
    constexpr auto operator[](ddc::DiscreteDomain<QueryDDims...> const& odomain)
    {
        return get_slice(odomain, std::make_index_sequence<base_type::NDims> {});
    }
};

template <
        class ElementType,
        class Domain,
        class NDTag,
        class LayoutStridedPolicy = std::experimental::layout_right,
        class MemorySpace = Kokkos::DefaultHostExecutionSpace::memory_space>
using VectorFieldView
        = VectorFieldSpan<const ElementType, Domain, NDTag, LayoutStridedPolicy, MemorySpace>;
