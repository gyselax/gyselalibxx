

# File multipatch\_field\_mem.hpp

[**File List**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_field\_mem.hpp**](multipatch__field__mem_8hpp.md)

[Go to the documentation of this file](multipatch__field__mem_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "multipatch_field.hpp"
#include "multipatch_type.hpp"

template <class T>
inline constexpr bool enable_multipatch_field_mem = false;

template <class T>
inline constexpr bool is_multipatch_field_mem_v
        = enable_multipatch_field_mem<std::remove_const_t<std::remove_reference_t<T>>>;

template <template <typename P> typename T, class... Patches>
class MultipatchFieldMem : public MultipatchType<T, Patches...>
{
    static_assert(
            (has_data_access_methods_v<T<Patches>> && ...),
            "The MultipatchFieldMem type should only contain instances of objects that can be "
            "manipulated like fields.");
    static_assert(
            (is_mem_type_v<T<Patches>> && ...),
            "The MultipatchFieldMem type should only contain instances of objects that allocate "
            "memory.");

public:
    using base_type = MultipatchType<T, Patches...>;

    using typename base_type::PatchOrdering;

    template <class Patch>
    using InternalIdxRangeOnPatch = typename T<Patch>::discrete_domain_type;

    template <class Patch>
    using InternalFieldOnPatch = typename T<Patch>::span_type;

    template <class Patch>
    using InternalConstFieldOnPatch = typename T<Patch>::view_type;

    template <template <typename P> typename OtherType, class... OPatches>
    friend class MultipatchFieldMem;

public:
    using span_type = MultipatchField<InternalFieldOnPatch, Patches...>;
    using view_type = MultipatchField<InternalConstFieldOnPatch, Patches...>;
    using discrete_domain_type = MultipatchType<InternalIdxRangeOnPatch, Patches...>;
    using memory_space = typename base_type::example_element::memory_space;
    using element_type = typename base_type::example_element::element_type;

public:
    explicit MultipatchFieldMem(T<Patches>... args) : base_type(args...) {}

    template <class MultipatchObj>
    explicit MultipatchFieldMem(MultipatchObj& other)
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    template <template <typename P> typename OtherType, class... OPatches>
    MultipatchFieldMem(MultipatchFieldMem<OtherType, OPatches...>&& other) : base_type(other)
    {
        static_assert(
                std::is_same_v<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<OPatches...>>,
                "Cannot create a MultipatchFieldMem from a temporary MultipatchFieldMem with a "
                "different "
                "ordering");
        static_assert(
                std::is_same_v<std::tuple<T<Patches>...>, std::tuple<OtherType<OPatches>...>>,
                "MultipatchFieldMems are not equivalent");
    }

    ~MultipatchFieldMem() noexcept = default;

    template <class Patch>
    auto get() const
    {
        return ::get_const_field(std::get<T<Patch>>(base_type::m_tuple));
    }

    template <class Patch>
    auto get()
    {
        return ::get_field(std::get<T<Patch>>(base_type::m_tuple));
    }

    auto idx_range() const
    {
        return MultipatchType<InternalIdxRangeOnPatch, Patches...>(
                get_idx_range(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    auto get_field()
    {
        return MultipatchField<InternalFieldOnPatch, Patches...>(
                ::get_field(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    auto span_view()
    {
        return get_field();
    }

    auto get_const_field() const
    {
        return MultipatchField<InternalConstFieldOnPatch, Patches...>(
                ::get_const_field(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    auto span_cview() const
    {
        return get_const_field();
    }
};

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_multipatch_type<MultipatchFieldMem<T, Patches...>> = true;

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_multipatch_field_mem<MultipatchFieldMem<T, Patches...>> = true;

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_mem_type<MultipatchFieldMem<T, Patches...>> = true;

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_data_access_methods<MultipatchFieldMem<T, Patches...>> = true;
```


