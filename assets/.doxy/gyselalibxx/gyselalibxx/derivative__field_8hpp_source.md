

# File derivative\_field.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**derivative\_field.hpp**](derivative__field_8hpp.md)

[Go to the documentation of this file](derivative__field_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "derivative_field_common.hpp"

template <class, class, int, class>
class DerivFieldMem;

template <
        class ElementType,
        class SupportType,
        class MemorySpace = Kokkos::HostSpace,
        class LayoutStridedPolicy = Kokkos::layout_right>
class DerivField;

template <class ElementType, class SupportType, class MemorySpace, class LayoutStridedPolicy>
inline constexpr bool enable_deriv_field<
        DerivField<ElementType, SupportType, MemorySpace, LayoutStridedPolicy>> = true;

template <class ElementType, class SupportType, class MemorySpace, class LayoutStridedPolicy>
inline constexpr bool enable_borrowed_deriv_field<
        DerivField<ElementType, SupportType, MemorySpace, LayoutStridedPolicy>> = true;

template <class ElementType, class SupportType, class MemorySpace, class LayoutStridedPolicy>
inline constexpr bool enable_data_access_methods<
        DerivField<ElementType, SupportType, MemorySpace, LayoutStridedPolicy>> = true;

namespace ddcHelper {

template <
        class FieldDst,
        class FieldSrc,
        std::enable_if_t<
                is_borrowed_deriv_field_v<FieldDst> && is_borrowed_deriv_field_v<FieldSrc>,
                bool> = true>
auto deepcopy(FieldDst&& dst, FieldSrc&& src)
{
    assert(dst.get_values_field().domain().extents() == src.get_values_field().domain().extents());

    DerivField dst_field = get_field(dst);
    DerivField src_field = get_field(src);

    dst_field.deepcopy(src_field);

    return dst_field;
}

template <class ExecSpace, class FieldDst, class FieldSrc>
auto deepcopy(ExecSpace const& execution_space, FieldDst&& dst, FieldSrc&& src)
{
    static_assert(is_borrowed_deriv_field_v<FieldDst>);
    static_assert(is_borrowed_deriv_field_v<FieldSrc>);

    assert(dst.get_values_field().idx_range().extents()
           == src.get_values_field().idx_range().extents());

    DerivField dst_field = get_field(dst);
    DerivField src_field = get_field(src);

    dst_field.deepcopy(execution_space, src_field);

    return dst_field;
}

} // namespace ddcHelper

template <class ElementType, class... DDims, class MemorySpace, class LayoutStridedPolicy>
class DerivField<ElementType, IdxRange<DDims...>, MemorySpace, LayoutStridedPolicy>
    : public DerivFieldCommon<
              Field<ElementType, IdxRange<DDims...>, MemorySpace, LayoutStridedPolicy>,
              IdxRange<DDims...>>
{
private:
    using base_type = DerivFieldCommon<
            Field<ElementType, IdxRange<DDims...>, MemorySpace, LayoutStridedPolicy>,
            IdxRange<DDims...>>;

public:
    using element_type = typename base_type::element_type;

    using discrete_domain_type = typename base_type::discrete_domain_type;
    using index_range_type = typename base_type::index_range_type;

    using discrete_element_type = typename base_type::discrete_element_type;

    using deriv_tags = typename base_type::deriv_tags;

    using physical_deriv_grids = typename detail::strip_deriv_t<deriv_tags>;

    using physical_grids = typename base_type::physical_grids;

    using physical_idx_range_type = typename base_type::physical_idx_range_type;

    using chunk_type = typename base_type::chunk_type;

    using span_type = DerivField<ElementType, IdxRange<DDims...>, MemorySpace, LayoutStridedPolicy>;

    using view_type
            = DerivField<ElementType const, IdxRange<DDims...>, MemorySpace, LayoutStridedPolicy>;

private:
    using discrete_deriv_idx_range_type = typename base_type::discrete_deriv_idx_range_type;

    using discrete_deriv_index_type = typename base_type::discrete_deriv_index_type;

    using reference = typename chunk_type::reference;

    using base_type::n_fields;

private:
    template <std::size_t ArrayIndex, class Tag>
    KOKKOS_FUNCTION constexpr Idx<Tag> get_chunk_subidx_range_1d_idx()
    {
        // The Tag describes a dimension for which a derivative is defined
        if constexpr (ddc::in_tags_v<Tag, deriv_tags>) {
            // If the Field at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, deriv_tags>)) {
                return Idx<Tag>(1);
            }
            // If the Field at this index doesn't contain derivatives of this dimension
            else {
                return Idx<Tag>(0);
            }
        } else {
            // Empty IdxRange to be discarded
            return Idx<Tag>();
        }
    }

    template <std::size_t ArrayIndex>
    KOKKOS_FUNCTION constexpr discrete_deriv_index_type get_chunk_subidx_range_idx()
    {
        return discrete_deriv_index_type(get_chunk_subidx_range_1d_idx<ArrayIndex, DDims>()...);
    }

    template <class Field, std::size_t... ArrayIndex>
    KOKKOS_FUNCTION constexpr void initialise_fields(
            Field const& chunks,
            std::index_sequence<ArrayIndex...>)
    {
        ((base_type::internal_fields[ArrayIndex] = chunks.internal_fields[chunks.get_array_index(
                  get_chunk_subidx_range_idx<ArrayIndex>())]),
         ...);
    }

    auto get_kokkos_view_from_internal_chunk(int index)
    {
        typename base_type::internal_mdspan_type field = base_type::internal_fields[index];
        auto kokkos_layout = ddc::detail::build_kokkos_layout(
                field.extents(),
                field.mapping(),
                std::make_index_sequence<sizeof...(DDims)> {});
        return Kokkos::View<
                ddc::detail::mdspan_to_kokkos_element_t<ElementType, sizeof...(DDims)>,
                decltype(kokkos_layout),
                MemorySpace>(field.data_handle(), kokkos_layout);
    }

public:
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField(DerivField const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField(DerivField&& other) = default;

    template <
            class OElementType,
            int NDerivs,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    explicit constexpr DerivField(
            DerivFieldMem<OElementType, index_range_type, NDerivs, Allocator>& field)
        : base_type(
                field.m_physical_idx_range,
                field.m_deriv_idx_range,
                field.m_cross_derivative_idx_range)
    {
        initialise_fields(field, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    // Disabled by SFINAE if `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            int NDerivs,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    explicit constexpr DerivField(
            DerivFieldMem<OElementType, index_range_type, NDerivs, Allocator> const& field)
        : base_type(
                field.m_physical_idx_range,
                field.m_deriv_idx_range,
                field.m_cross_derivative_idx_range)
    {
        initialise_fields(field, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    template <class OElementType>
    KOKKOS_FUNCTION constexpr DerivField(
            DerivField<OElementType, index_range_type, MemorySpace, LayoutStridedPolicy> const&
                    field)
        : base_type(
                field.m_physical_idx_range,
                field.m_deriv_idx_range,
                field.m_cross_derivative_idx_range)
    {
        initialise_fields(field, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    KOKKOS_DEFAULTED_FUNCTION ~DerivField() = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField& operator=(DerivField const& other) = default;

    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField& operator=(DerivField&& other) = default;

    template <class OElementType, class OLayoutStridedPolicy, class OMemorySpace>
    void deepcopy(
            DerivField<OElementType, index_range_type, OMemorySpace, OLayoutStridedPolicy> src)
    {
        for (int i(0); i < n_fields; ++i) {
            auto kokkos_span = get_kokkos_view_from_internal_chunk(i);
            auto src_kokkos_span = src.get_kokkos_view_from_internal_chunk(i);
            Kokkos::deep_copy(kokkos_span, src_kokkos_span);
        }
    }

    template <class ExecSpace, class OElementType, class OMemorySpace, class OLayoutStridedPolicy>
    void deepcopy(
            ExecSpace const& execution_space,
            DerivField<OElementType, index_range_type, OMemorySpace, OLayoutStridedPolicy> src)
    {
        for (int i(0); i < n_fields; ++i) {
            auto kokkos_span = get_kokkos_view_from_internal_chunk(i);
            auto src_kokkos_span = src.get_kokkos_view_from_internal_chunk(i);
            Kokkos::deep_copy(execution_space, kokkos_span, src_kokkos_span);
        }
    }

    template <class... DElem>
    KOKKOS_FUNCTION constexpr reference operator()(DElem... elems) const noexcept
    {
        static_assert((ddc::is_discrete_element_v<DElem> && ...));
        using full_index_type = detail::combine_t<DElem...>;
        full_index_type elem(elems...);
        return base_type::get_internal_field(elem)();
    }

    KOKKOS_FUNCTION constexpr view_type span_cview() const
    {
        return view_type(*this);
    }

    KOKKOS_FUNCTION constexpr span_type span_view() const
    {
        return *this;
    }
};

template <
        class ElementType,
        class SupportType,
        class MemorySpace = Kokkos::HostSpace,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DerivConstField
        = DerivField<ElementType const, SupportType, MemorySpace, LayoutStridedPolicy>;

namespace detail {
template <
        class NewMemorySpace,
        class ElementType,
        class SupportType,
        class MemorySpace,
        class Layout>
struct OnMemorySpace<NewMemorySpace, DerivField<ElementType, SupportType, MemorySpace, Layout>>
{
    using type = DerivField<ElementType, SupportType, NewMemorySpace, Layout>;
};
}; // namespace detail
```


