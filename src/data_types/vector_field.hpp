// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"


template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
class VectorField;

template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace,
        class LayoutStridedPolicy>
inline constexpr bool enable_vector_field<VectorField<
        ElementType,
        IdxRangeType,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>> = true;

template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace,
        class LayoutStridedPolicy>
inline constexpr bool enable_data_access_methods<VectorField<
        ElementType,
        IdxRangeType,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>> = true;

template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace,
        class LayoutStridedPolicy>
inline constexpr bool enable_borrowed_vector_field<VectorField<
        ElementType,
        IdxRangeType,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>> = true;


/**
 * @brief A class which holds multiple (scalar) fields in order to represent a vector field.
 *
 * @tparam ElementType The data type of a scalar element of the vector field.
 * @tparam IdxRangeType
 * @tparam VectorIndexSetType A VectorIndexSet describing the dimensions described by the scalar elements of a vector field element.
 * @tparam MemorySpace The memory space (CPU/GPU).
 * @tparam LayoutStridedPolicy The memory layout. See DDC.
 */
template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace,
        class LayoutStridedPolicy>
class VectorField
    : public VectorFieldCommon<
              Field<ElementType, IdxRangeType, MemorySpace, LayoutStridedPolicy>,
              VectorIndexSetType>
{
public:
    /**
     * @brief Type describing the object which can be extracted from this VectorField using the get<> function.
     */
    using field_type = Field<ElementType, IdxRangeType, MemorySpace, LayoutStridedPolicy>;

private:
    using base_type = VectorFieldCommon<field_type, VectorIndexSetType>;

public:
    /// The type of an element in one of the Fields comprising the VectorField
    using element_type = typename base_type::element_type;

    using typename base_type::NDTypeTag;

    using typename base_type::chunk_span_type;
    using typename base_type::chunk_view_type;

public:
    /**
     * @brief A type which can hold a modifiable reference to a VectorFieldMem.
     *
     * A type which can hold a reference to a VectorFieldMem. If this object is modifiable then
     * so is the span type.
     * This is a DDC keyword used to make this class interchangeable with Field.
     */
    using span_type = VectorField<
            ElementType,
            IdxRangeType,
            VectorIndexSetType,
            MemorySpace,
            LayoutStridedPolicy>;
    /**
     * @brief A type which can hold a constant reference to a VectorFieldMem.
     * This is a DDC keyword used to make this class interchangeable with Field.
     */
    using view_type = VectorField<
            const ElementType,
            IdxRangeType,
            VectorIndexSetType,
            MemorySpace,
            LayoutStridedPolicy>;

    /**
     * @brief Type describing the way in which the data is laid out in the Field memory.
     *
     * Type describing the way in which the data is laid out in the Field memory.
     * I.e. it describes whether it is contiguous or not.
     */
    using layout_type = LayoutStridedPolicy;

    /**
     * @brief The type of the index range on which the field is defined.
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC IdxRange types are referred to as DiscreteDomain types.
     */
    using discrete_domain_type = typename base_type::discrete_domain_type;
    /// @brief The IdxRange on which the fields in this object are defined.
    using index_range_type = discrete_domain_type;

    /**
     * @brief The type of the memory space where the field is saved (CPU vs GPU).
     */
    using memory_space = typename field_type::memory_space;

private:
    template <class, class, class, class, class>
    friend class VectorField;

    template <class OElementType, class Allocator, std::size_t... Is>
    KOKKOS_FUNCTION constexpr VectorField(
            VectorFieldMem<OElementType, IdxRangeType, VectorIndexSetType, Allocator>& other,
            std::index_sequence<Is...> const&) noexcept
        : base_type((field_type(
                ddcHelper::get<ddc::type_seq_element_t<Is, VectorIndexSetType>>(other)))...)
    {
    }

    /** Constructs a new VectorField from a VectorFieldMem, yields a new view to the same data
     * @param other the VectorFieldMem to view
     */
    // Disabled by SFINAE in the case of `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            std::size_t... Is,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            class Allocator>
    KOKKOS_FUNCTION constexpr VectorField(
            VectorFieldMem<OElementType, IdxRangeType, VectorIndexSetType, Allocator> const& other,
            std::index_sequence<Is...> const&) noexcept
        : base_type((field_type(
                ddcHelper::get<ddc::type_seq_element_t<Is, VectorIndexSetType>>(other)))...)
    {
    }

    /** Constructs a new VectorField by copy of a chunk, yields a new view to the same data
     * @param other the VectorField to move
     */
    template <class OElementType, std::size_t... Is>
    KOKKOS_FUNCTION constexpr VectorField(
            VectorField<
                    OElementType,
                    index_range_type,
                    VectorIndexSetType,
                    MemorySpace,
                    LayoutStridedPolicy> const& other,
            std::index_sequence<Is...> const&) noexcept
        : base_type((field_type(
                ddcHelper::get<ddc::type_seq_element_t<Is, VectorIndexSetType>>(other)))...)
    {
    }

    template <class SliceType, std::size_t... Is>
    constexpr auto get_slice(SliceType const& slice_spec, std::index_sequence<Is...> const&)
    {
        auto chunk_slices = std::make_tuple(
                this->template get<
                        ddc::type_seq_element_t<Is, VectorIndexSetType>>()[slice_spec]...);
        using FieldType = std::tuple_element_t<0, decltype(chunk_slices)>;
        return VectorField<
                ElementType,
                typename FieldType::discrete_domain_type,
                VectorIndexSetType,
                typename FieldType::memory_space,
                typename FieldType::layout_type>(std::move(std::get<Is>(chunk_slices))...);
    }

    /** Element access using a multi-dimensional Idx
     * @param delems discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims, typename T, T... ints>
    KOKKOS_FUNCTION element_type const operator()(
            Idx<ODDims...> const& delems,
            std::integer_sequence<T, ints...>) const noexcept
    {
        return element_type((base_type::m_values[ints](delems))...);
    }

public:
    /// Empty VectorField
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField() = default;

    /// VectorField destructor
    KOKKOS_DEFAULTED_FUNCTION ~VectorField() = default;

    /** Constructs a new VectorField by copy, yields a new view to the same data
     * @param other the VectorField to copy
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField(VectorField const& other) = default;

    /** Constructs a new VectorField by move
     * @param other the VectorField to move
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField(VectorField&& other) = default;

    /** Constructs a new VectorField from a VectorFieldMem, yields a new view to the same data
     * @param other the VectorFieldMem to view
     */
    template <class OElementType, class Allocator>
    KOKKOS_FUNCTION explicit constexpr VectorField(
            VectorFieldMem<OElementType, IdxRangeType, VectorIndexSetType, Allocator>&
                    other) noexcept
        : VectorField(other, std::make_index_sequence<base_type::NDims> {})
    {
    }

    /** Constructs a new VectorField from a VectorFieldMem, yields a new view to the same data
     * @param other the VectorFieldMem to view
     */
    // Disabled by SFINAE in the case of `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            class Allocator>
    KOKKOS_FUNCTION explicit constexpr VectorField(
            VectorFieldMem<OElementType, IdxRangeType, VectorIndexSetType, Allocator> const&
                    other) noexcept
        : VectorField(other, std::make_index_sequence<base_type::NDims> {})
    {
    }

    template <class OElementType, class Allocator>
    VectorField(VectorFieldMem<OElementType, IdxRangeType, VectorIndexSetType, Allocator>&& other)
            = delete;

    /** Constructs a new VectorField by copy of a chunk, yields a new view to the same data
     * @param other the VectorField to move
     */
    template <class OElementType>
    KOKKOS_FUNCTION constexpr VectorField(VectorField<
                                          OElementType,
                                          index_range_type,
                                          VectorIndexSetType,
                                          MemorySpace,
                                          LayoutStridedPolicy> const& other) noexcept
        : VectorField(other, std::make_index_sequence<base_type::NDims> {})
    {
    }

    /** Constructs a new VectorField from scratch
     * @param ptr the allocation pointer to the data
     * @param idx_range the index range that sustains the view
     */
    template <
            class... OElementType,
            class = std::enable_if_t<
                    std::conjunction_v<std::is_same<OElementType, ElementType>...>>,
            class = std::enable_if_t<sizeof...(OElementType) == base_type::NDims>>
    KOKKOS_FUNCTION VectorField(index_range_type const& idx_range, OElementType*... ptr)
        : base_type((field_type(ptr, idx_range))...)
    {
    }

    /** Constructs a new VectorField containing references to Field.
     * @param fields The Fields.
     */
    template <
            class... FieldType,
            class = std::enable_if_t<std::conjunction_v<std::is_same<FieldType, field_type>...>>>
    KOKKOS_FUNCTION constexpr VectorField(FieldType... fields) : base_type(std::move(fields)...)
    {
    }

    /** Copy-assigns a new value to this VectorField, yields a new view to the same data
     * @param other the VectorField to copy
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField& operator=(VectorField const& other) = default;

    /** Move-assigns a new value to this VectorField
     * @param other the VectorField to move
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField& operator=(VectorField&& other) = default;

    /**
     * Get a constant reference to the vector field referred to by this vector field span.
     *
     * This function is designed to match the equivalent function in DDC. In Gysela it should
     * not be called directly. Instead the global function get_const_field should be used.
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
     * This function is designed to match the equivalent function in DDC. In Gysela it should
     * not be called directly. Instead the global function get_field should be used.
     *
     * @return A constant reference to the vector field.
     */
    constexpr span_type span_view() const
    {
        return *this;
    }

    /** Element access using a list of Idxs
     * @param delems 1D discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims>
    KOKKOS_FUNCTION element_type const operator()(
            ddc::DiscreteElement<ODDims> const&... delems) const noexcept
    {
        Idx<ODDims...> delem_idx(delems...);
        return this->
        operator()(delem_idx, std::make_integer_sequence<int, element_type::size()> {});
    }

    /** Element access using a multi-dimensional Idx
     * @param delems discrete coordinates
     * @return copy of this element
     */
    template <class... ODDims, class = std::enable_if_t<sizeof...(ODDims) != 1>>
    KOKKOS_FUNCTION element_type const operator()(Idx<ODDims...> const& delems) const noexcept
    {
        return this->operator()(delems, std::make_integer_sequence<int, element_type::size()> {});
    }


    /**
     * @brief Slice out some dimensions.
     *
     * Get the VectorFieldMem on the reduced index range which is obtained by indexing
     * the dimensions QueryDDims at the position slice_spec.
     *
     * @param[in] slice_spec The slice describing the index range of interest.
     *
     * @return A reference to the vector field on the sliced index range.
     */
    template <class... QueryDDims>
    constexpr auto operator[](Idx<QueryDDims...> const& slice_spec)
    {
        return get_slice(slice_spec, std::make_index_sequence<base_type::NDims> {});
    }

    /**
     * @brief Slice out some dimensions.
     *
     * Get the VectorFieldMem on the reduced index range passed as an argument.
     *
     * @param[in] oidx_range The index range of interest.
     *
     * @return A reference to the vector field on the sliced index range.
     */
    template <class... QueryDDims>
    constexpr auto operator[](IdxRange<QueryDDims...> const& oidx_range)
    {
        return get_slice(oidx_range, std::make_index_sequence<base_type::NDims> {});
    }

    /**
     * @brief Get the Field describing the component in the QueryTag direction.
     *
     * @return The field in the specified direction.
     */
    template <class QueryTag>
    inline constexpr chunk_span_type get() const noexcept
    {
        static_assert(ddc::in_tags_v<QueryTag, NDTypeTag>, "requested Tag absent from Vector");
        return base_type::m_values[ddc::type_seq_rank_v<QueryTag, NDTypeTag>].span_view();
    }
};

template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using VectorConstField = VectorField<
        const ElementType,
        IdxRangeType,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>;

template <
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DVectorField
        = VectorField<double, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy>;

template <
        class IdxRangeType,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DVectorConstField = VectorConstField<
        double,
        IdxRangeType,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>;

namespace detail {

/**
 * @brief Get a new `VectorField` type with the same parametrisation
 * except in the memory space which is set to NewMemorySpace.
 * @tparam NewMemorySpace The new memory space. 
 * @tparam ElementType Type of the elements in the ddc::Chunk of the VectorFieldMem.
 * @tparam SupportType Type of the domain of the ddc::Chunk in the VectorFieldMem.
 * @tparam VectorIndexSetType VectorIndexSet object storing directions of the VectorFieldMem as dimensions. 
 *               The dimensions refer to the dimensions of the arrival domain of the VectorFieldMem. 
 * @tparam Layout Layout tag (see Kokkos).
 * @tparam MemorySpace The original memory space of the chunk of the VectorFieldMem.
 * @see VectorField
 */
template <
        class NewMemorySpace,
        class ElementType,
        class SupportType,
        class VectorIndexSetType,
        class MemorySpace,
        class Layout>
struct OnMemorySpace<
        NewMemorySpace,
        VectorField<ElementType, SupportType, VectorIndexSetType, MemorySpace, Layout>>
{
    using type = VectorField<ElementType, SupportType, VectorIndexSetType, NewMemorySpace, Layout>;
};

} // namespace detail

namespace ddcHelper {

template <
        class ExecSpace,
        class ElementType,
        class IdxRangeType,
        class... Dims,
        class MemorySpace,
        class LayoutStridedPolicy>
auto create_mirror_view_and_copy(
        ExecSpace exec_space,
        VectorField<
                ElementType,
                IdxRangeType,
                VectorIndexSet<Dims...>,
                MemorySpace,
                LayoutStridedPolicy> field)
{
    if constexpr (Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible) {
        return field;
    } else {
        VectorFieldMem<
                std::remove_const_t<ElementType>,
                IdxRangeType,
                VectorIndexSet<Dims...>,
                typename ExecSpace::memory_space>
                field_alloc(get_idx_range(field));
        ((ddc::parallel_deepcopy(field_alloc.template get<Dims>(), field.template get<Dims>())),
         ...);
        return field_alloc;
    }
}

/**
 * @brief Copy the elements of a vector into a vector field at a given index.
 *
 * @tparam ElementType Type of the elements in the ddc::Chunk of the VectorField and Vector.
 * @tparam IdxRangeType Type of the IndexRange over which the VectorField is defined. 
 * @tparam Dims... Dimensions of the VectorField and Vector. 
 * @tparam MemorySpace The original memory space of the chunk of the VectorFieldMem.
 * @tparam LayoutStridedPolicy Layout tag (see Kokkos).
 *
 * @param[out] field On output, a VectorField containing the values of vector 
 * at the index idx.
 * @param[in] idx An index so specify where the values of the Vector should be 
 * copied in the VectorField. 
 * @param[in] vector A Vector to be copied in the VectorField at the index idx.
 */
template <
        class ElementType,
        class IdxRangeType,
        class... Dims,
        class MemorySpace,
        class LayoutStridedPolicy>
KOKKOS_INLINE_FUNCTION void assign_vector_field_element(
        VectorField<
                ElementType,
                IdxRangeType,
                VectorIndexSet<Dims...>,
                MemorySpace,
                LayoutStridedPolicy> field,
        typename IdxRangeType::discrete_element_type idx,
        Vector<ElementType, Dims...> vector)
{
    ((ddcHelper::get<Dims>(field)(idx) = ddcHelper::get<Dims>(vector)), ...);
}
} // namespace ddcHelper
