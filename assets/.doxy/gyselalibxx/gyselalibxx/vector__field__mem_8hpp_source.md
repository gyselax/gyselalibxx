

# File vector\_field\_mem.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**vector\_field\_mem.hpp**](vector__field__mem_8hpp.md)

[Go to the documentation of this file](vector__field__mem_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "vector_field_common.hpp"
#include "vector_index_tools.hpp"

template <
        class ElementType,
        class IdxRangeType,
        class VectorIndexSetType,
        class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
class VectorFieldMem;


template <class ElementType, class IdxRangeType, class DimSeq, class MemSpace>
inline constexpr bool
        enable_vector_field<VectorFieldMem<ElementType, IdxRangeType, DimSeq, MemSpace>> = true;

template <class ElementType, class IdxRangeType, class DimSeq, class Allocator>
inline constexpr bool enable_data_access_methods<
        VectorFieldMem<ElementType, IdxRangeType, DimSeq, Allocator>> = true;

template <class ElementType, class IdxRangeType, class DimSeq, class Allocator>
inline constexpr bool
        enable_mem_type<VectorFieldMem<ElementType, IdxRangeType, DimSeq, Allocator>> = true;

template <class, class, class, class, class>
class VectorField;

template <class ElementType, class IdxRangeType, class VectorIndexSetType, class MemSpace>
class VectorFieldMem
    : public VectorFieldCommon<FieldMem<ElementType, IdxRangeType, MemSpace>, VectorIndexSetType>
{
public:
    using chunk_type = FieldMem<ElementType, IdxRangeType, MemSpace>;

private:
    using base_type = VectorFieldCommon<chunk_type, VectorIndexSetType>;

public:
    using typename base_type::element_type;

    using typename base_type::NDTypeTag;

    using typename base_type::chunk_span_type;
    using typename base_type::chunk_view_type;

    using Allocator = ddc::KokkosAllocator<ElementType, MemSpace>;

public:
    using span_type = VectorField<
            ElementType,
            IdxRangeType,
            VectorIndexSetType,
            MemSpace,
            Kokkos::layout_right>;

    using view_type = VectorField<
            const ElementType,
            IdxRangeType,
            VectorIndexSetType,
            MemSpace,
            Kokkos::layout_right>;

    using discrete_domain_type = typename base_type::discrete_domain_type;
    using index_range_type = discrete_domain_type;

    using memory_space = typename chunk_type::memory_space;

private:
    template <std::size_t... Is>
    explicit VectorFieldMem(
            index_range_type const& idx_range,
            Allocator allocator,
            std::index_sequence<Is...> const&)
        : base_type(((void)Is, chunk_type(idx_range, allocator))...)
    {
    }

    template <class... ODDims, typename T, T... ints>
    element_type const operator()(Idx<ODDims...> const& delems, std::integer_sequence<T, ints...>)
            const noexcept
    {
        return element_type((base_type::m_values[ints](delems))...);
    }

public:
    VectorFieldMem() = default;

    explicit VectorFieldMem(index_range_type const& idx_range, Allocator allocator = Allocator())
        : VectorFieldMem(idx_range, allocator, std::make_index_sequence<base_type::NDims> {})
    {
    }

    VectorFieldMem(VectorFieldMem const& other) = delete;

    VectorFieldMem(VectorFieldMem&& other) = default;

    VectorFieldMem& operator=(VectorFieldMem const& other) = delete;

    VectorFieldMem& operator=(VectorFieldMem&& other) = default;

    ~VectorFieldMem() = default;

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

    template <class... ODDims>
    element_type const operator()(ddc::DiscreteElement<ODDims> const&... delems) const noexcept
    {
        Idx<ODDims...> delem_idx(delems...);
        return this->
        operator()(delem_idx, std::make_integer_sequence<int, element_type::size()> {});
    }

    template <class... ODDims, class = std::enable_if_t<sizeof...(ODDims) != 1>>
    element_type const operator()(Idx<ODDims...> const& delems) const noexcept
    {
        return this->operator()(delems, std::make_integer_sequence<int, element_type::size()> {});
    }


    template <class... QueryDDims>
    auto operator[](Idx<QueryDDims...> const& slice_spec) const
    {
        return span_cview()[slice_spec];
    }

    template <class... QueryDDims>
    auto operator[](Idx<QueryDDims...> const& slice_spec)
    {
        return span_view()[slice_spec];
    }

    template <class... QueryDDims>
    auto operator[](IdxRange<QueryDDims...> const& oidx_range) const
    {
        return span_cview()[oidx_range];
    }

    template <class... QueryDDims>
    auto operator[](IdxRange<QueryDDims...> const& oidx_range)
    {
        return span_view()[oidx_range];
    }

    template <class QueryTag>
    inline constexpr chunk_span_type get() noexcept
    {
        static_assert(ddc::in_tags_v<QueryTag, NDTypeTag>, "requested Tag absent from Vector");
        return base_type::m_values[ddc::type_seq_rank_v<QueryTag, NDTypeTag>].span_view();
    }

    template <class QueryTag>
    inline constexpr chunk_view_type get() const noexcept
    {
        static_assert(ddc::in_tags_v<QueryTag, NDTypeTag>, "requested Tag absent from Vector");
        return base_type::m_values[ddc::type_seq_rank_v<QueryTag, NDTypeTag>].span_cview();
    }
};

template <
        class IdxRangeType,
        class VectorIndexSetType,
        class MemSpace = Kokkos::DefaultExecutionSpace::memory_space>
using DVectorFieldMem = VectorFieldMem<double, IdxRangeType, VectorIndexSetType, MemSpace>;

namespace detail {
template <
        class NewMemorySpace,
        class ElementType,
        class SupportType,
        class VectorIndexSetType,
        class MemSpace>
struct OnMemorySpace<
        NewMemorySpace,
        VectorFieldMem<ElementType, SupportType, VectorIndexSetType, MemSpace>>
{
    using type = VectorFieldMem<ElementType, SupportType, VectorIndexSetType, NewMemorySpace>;
};
} // namespace detail
```


