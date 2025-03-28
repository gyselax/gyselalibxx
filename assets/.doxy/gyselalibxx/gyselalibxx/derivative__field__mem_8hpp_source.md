

# File derivative\_field\_mem.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**derivative\_field\_mem.hpp**](derivative__field__mem_8hpp.md)

[Go to the documentation of this file](derivative__field__mem_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "derivative_field.hpp"
#include "derivative_field_common.hpp"

template <class ElementType, class Domain, int NDerivs, class MemSpace = Kokkos::HostSpace>
class DerivFieldMem;

template <class ElementType, class SupportType, int NDerivs, class MemSpace>
inline constexpr bool
        enable_deriv_field<DerivFieldMem<ElementType, SupportType, NDerivs, MemSpace>> = true;

template <class ElementType, class SupportType, int NDerivs, class Allocator>
inline constexpr bool enable_data_access_methods<
        DerivFieldMem<ElementType, SupportType, NDerivs, Allocator>> = true;

template <class ElementType, class SupportType, int NDerivs, class Allocator>
inline constexpr bool
        enable_mem_type<DerivFieldMem<ElementType, SupportType, NDerivs, Allocator>> = true;


template <class ElementType, class... DDims, int NDerivs, class MemSpace>
class DerivFieldMem<ElementType, IdxRange<DDims...>, NDerivs, MemSpace>
    : public DerivFieldCommon<
              FieldMem<ElementType, IdxRange<DDims...>, MemSpace>,
              IdxRange<DDims...>>
{
private:
    using base_type = DerivFieldCommon<
            FieldMem<ElementType, IdxRange<DDims...>, MemSpace>,
            IdxRange<DDims...>>;

public:
    using element_type = typename base_type::element_type;

    using discrete_domain_type = typename base_type::discrete_domain_type;

    using discrete_element_type = typename base_type::discrete_element_type;

    using deriv_tags = typename base_type::deriv_tags;

    using physical_deriv_grids = typename base_type::physical_deriv_grids;

    using physical_grids = typename base_type::physical_grids;

    using physical_idx_range_type = typename base_type::physical_idx_range_type;

    using chunk_type = typename base_type::chunk_type;

    template <class, class, class, class>
    friend class DerivField;

    using span_type = DerivField<
            ElementType,
            IdxRange<DDims...>,
            typename chunk_type::memory_space,
            typename chunk_type::layout_type>;

    using view_type = DerivField<
            ElementType const,
            IdxRange<DDims...>,
            typename chunk_type::memory_space,
            typename chunk_type::layout_type>;

private:
    using physical_deriv_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<physical_deriv_grids>;

    using discrete_deriv_idx_range_type = typename base_type::discrete_deriv_idx_range_type;

    using discrete_deriv_index_type = typename base_type::discrete_deriv_index_type;

    using mapping_type = typename chunk_type::mapping_type;

    using extents_type = typename chunk_type::extents_type;

    using allocator_type = ddc::KokkosAllocator<element_type, typename chunk_type::memory_space>;

    using internal_mdspan_type = typename base_type::internal_mdspan_type;

    using base_type::n_fields;

private:
    template <std::size_t ArrayIndex, class Tag>
    IdxRange<Tag> get_idx_range(
            physical_idx_range_type val_idx_range,
            physical_deriv_idx_range_type deriv_idx_ranges)
    {
        // If the Tag describes a dimension for which a derivative is defined
        if constexpr (ddc::in_tags_v<Tag, physical_deriv_grids>) {
            // If the FieldMem at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, physical_deriv_grids>)) {
                return ddc::select<Tag>(deriv_idx_ranges);
            }
            // If the FieldMem at this index doesn't contain derivatives of this dimension
            else {
                return ddc::select<Tag>(val_idx_range);
            }
        }
        // If the Tag describes a derivative
        else if constexpr (ddc::in_tags_v<Tag, deriv_tags>) {
            // If the FieldMem at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, deriv_tags>)) {
                return IdxRange<Tag>(Idx<Tag>(1), IdxStep<Tag>(NDerivs));
            }
            // If the FieldMem at this index doesn't contain derivatives of this dimension
            else {
                return IdxRange<Tag>(Idx<Tag>(0), IdxStep<Tag>(1));
            }
        }
        // If the Tag describes a dimension for which derivatives are not defined.
        else {
            return ddc::select<Tag>(val_idx_range);
        }
    }

    template <std::size_t ArrayIndex>
    IdxRange<DDims...> get_chunk_idx_range(
            physical_idx_range_type val_idx_range,
            physical_deriv_idx_range_type deriv_idx_ranges)
    {
        return IdxRange<DDims...>(
                get_idx_range<ArrayIndex, DDims>(val_idx_range, deriv_idx_ranges)...);
    }

    template <class QueryDDim, std::size_t ArrayIndex>
    std::size_t get_mdspan_size()
    {
        if constexpr (ddc::in_tags_v<QueryDDim, deriv_tags>) {
            if (ArrayIndex & (1 << ddc::type_seq_rank_v<QueryDDim, deriv_tags>)) {
                return NDerivs;
            } else {
                return 1;
            }
        } else if constexpr (ddc::in_tags_v<QueryDDim, physical_deriv_grids>) {
            if constexpr (
                    ArrayIndex & (1 << ddc::type_seq_rank_v<ddc::Deriv<QueryDDim>, deriv_tags>)) {
                IdxRangeSlice<QueryDDim> idx_range_local(base_type::m_cross_derivative_idx_range);
                return idx_range_local.extents().value();
            } else {
                IdxRange<QueryDDim> idx_range_local(base_type::m_physical_idx_range);
                return idx_range_local.extents().value();
            }
        } else {
            IdxRange<QueryDDim> idx_range_local(base_type::m_physical_idx_range);
            return idx_range_local.extents().value();
        }
    }

    template <std::size_t ArrayIndex>
    std::enable_if_t<std::is_constructible_v<mapping_type, extents_type>, internal_mdspan_type>
    make_internal_mdspan(allocator_type allocator)
    {
        std::size_t alloc_size(((get_mdspan_size<DDims, ArrayIndex>()) * ...));
        element_type* ptr = allocator.allocate("", alloc_size);

        extents_type extents_r(get_mdspan_size<DDims, ArrayIndex>()...);
        mapping_type mapping_r(extents_r);

        std::array<std::size_t, sizeof...(DDims)> strides_s {
                mapping_r.stride(ddc::type_seq_rank_v<DDims, ddc::detail::TypeSeq<DDims...>>)...};
        Kokkos::layout_stride::mapping<extents_type> mapping_s(extents_r, strides_s);
        return internal_mdspan_type(ptr, mapping_s);
    }

    template <std::size_t... ArrayIndex>
    void initialise_chunks(allocator_type allocator, std::index_sequence<ArrayIndex...>)
    {
        ((base_type::internal_fields[ArrayIndex] = make_internal_mdspan<ArrayIndex>(allocator)),
         ...);
    }

public:
    template <class... DerivDoms>
    DerivFieldMem(
            physical_idx_range_type val_idx_range,
            IdxRangeSlice<DerivDoms>... m_deriv_idx_range)
        : DerivFieldMem(allocator_type(), val_idx_range, m_deriv_idx_range...)
    {
    }

    template <class... DerivDoms>
    DerivFieldMem(
            allocator_type allocator,
            physical_idx_range_type val_idx_range,
            IdxRangeSlice<DerivDoms>... m_deriv_idx_range)
        : base_type(
                val_idx_range,
                discrete_deriv_idx_range_type(IdxRange<ddc::Deriv<DerivDoms>>(
                        Idx<ddc::Deriv<DerivDoms>>(1),
                        IdxStep<ddc::Deriv<DerivDoms>>(NDerivs))...),
                to_subidx_range_collection<physical_deriv_grids>(m_deriv_idx_range...))
    {
        static_assert(
                ddc::type_seq_same_v<ddc::detail::TypeSeq<DerivDoms...>, physical_deriv_grids>);
        initialise_chunks(allocator, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    ~DerivFieldMem() = default;

    DerivFieldMem& operator=(DerivFieldMem const& other) = delete;

    DerivFieldMem& operator=(DerivFieldMem&& other) = default;

    template <class... DElem>
    element_type& operator()(DElem... elems) noexcept
    {
        static_assert((ddc::is_discrete_element_v<DElem> && ...));
        using full_index_type = detail::combine_t<DElem...>;
        full_index_type elem(elems...);
        return base_type::get_internal_field(elem)();
    }

    template <class... DElem>
    element_type const& operator()(DElem... elems) const noexcept
    {
        static_assert((ddc::is_discrete_element_v<DElem> && ...));
        using full_index_type = detail::combine_t<DElem...>;
        full_index_type elem(elems...);
        return base_type::get_internal_field(elem)();
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

namespace detail {
template <class NewMemorySpace, class ElementType, class SupportType, int NDerivs, class MemSpace>
struct OnMemorySpace<NewMemorySpace, DerivFieldMem<ElementType, SupportType, NDerivs, MemSpace>>
{
    using type = DerivFieldMem<ElementType, SupportType, NDerivs, NewMemorySpace>;
};
} // namespace detail
```


