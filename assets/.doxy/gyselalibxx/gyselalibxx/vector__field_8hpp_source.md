

# File vector\_field.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**vector\_field.hpp**](vector__field_8hpp.md)

[Go to the documentation of this file](vector__field_8hpp.md)


```C++
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
    using field_type = Field<ElementType, IdxRangeType, MemorySpace, LayoutStridedPolicy>;

private:
    using base_type = VectorFieldCommon<field_type, VectorIndexSetType>;

public:
    using element_type = typename base_type::element_type;

    using typename base_type::NDTypeTag;

    using typename base_type::chunk_span_type;
    using typename base_type::chunk_view_type;

public:
    using span_type = VectorField<
            ElementType,
            IdxRangeType,
            VectorIndexSetType,
            MemorySpace,
            LayoutStridedPolicy>;
    using view_type = VectorField<
            const ElementType,
            IdxRangeType,
            VectorIndexSetType,
            MemorySpace,
            LayoutStridedPolicy>;

    using layout_type = LayoutStridedPolicy;

    using discrete_domain_type = typename base_type::discrete_domain_type;
    using index_range_type = discrete_domain_type;

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

    template <class... ODDims, typename T, T... ints>
    KOKKOS_FUNCTION element_type const operator()(
            Idx<ODDims...> const& delems,
            std::integer_sequence<T, ints...>) const noexcept
    {
        return element_type((base_type::m_values[ints](delems))...);
    }

public:
    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField() = default;

    KOKKOS_DEFAULTED_FUNCTION ~VectorField() = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField(VectorField const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField(VectorField&& other) = default;

    template <class OElementType, class Allocator>
    KOKKOS_FUNCTION explicit constexpr VectorField(
            VectorFieldMem<OElementType, IdxRangeType, VectorIndexSetType, Allocator>&
                    other) noexcept
        : VectorField(other, std::make_index_sequence<base_type::NDims> {})
    {
    }

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

    template <
            class... OElementType,
            class = std::enable_if_t<
                    std::conjunction_v<std::is_same<OElementType, ElementType>...>>,
            class = std::enable_if_t<sizeof...(OElementType) == base_type::NDims>>
    KOKKOS_FUNCTION VectorField(index_range_type const& idx_range, OElementType*... ptr)
        : base_type((field_type(ptr, idx_range))...)
    {
    }

    template <
            class... FieldType,
            class = std::enable_if_t<std::conjunction_v<std::is_same<FieldType, field_type>...>>>
    KOKKOS_FUNCTION constexpr VectorField(FieldType... fields) : base_type(std::move(fields)...)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField& operator=(VectorField const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr VectorField& operator=(VectorField&& other) = default;

    constexpr view_type span_cview() const
    {
        return view_type(*this);
    }

    constexpr span_type span_view() const
    {
        return *this;
    }

    template <class... ODDims>
    KOKKOS_FUNCTION element_type const operator()(
            ddc::DiscreteElement<ODDims> const&... delems) const noexcept
    {
        Idx<ODDims...> delem_idx(delems...);
        return this->
        operator()(delem_idx, std::make_integer_sequence<int, element_type::size()> {});
    }

    template <class... ODDims, class = std::enable_if_t<sizeof...(ODDims) != 1>>
    KOKKOS_FUNCTION element_type const operator()(Idx<ODDims...> const& delems) const noexcept
    {
        return this->operator()(delems, std::make_integer_sequence<int, element_type::size()> {});
    }


    template <class... QueryDDims>
    constexpr auto operator[](Idx<QueryDDims...> const& slice_spec)
    {
        return get_slice(slice_spec, std::make_index_sequence<base_type::NDims> {});
    }

    template <class... QueryDDims>
    constexpr auto operator[](IdxRange<QueryDDims...> const& oidx_range)
    {
        return get_slice(oidx_range, std::make_index_sequence<base_type::NDims> {});
    }

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
```


