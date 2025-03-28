

# File multipatch\_field.hpp

[**File List**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_field.hpp**](multipatch__field_8hpp.md)

[Go to the documentation of this file](multipatch__field_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include "multipatch_type.hpp"


template <class T>
inline constexpr bool enable_multipatch_field = false;

template <class T>
inline constexpr bool is_multipatch_field_v
        = enable_multipatch_field<std::remove_const_t<std::remove_reference_t<T>>>;

template <template <typename P> typename T, class... Patches>
class MultipatchField : public MultipatchType<T, Patches...>
{
    static_assert(
            (has_data_access_methods_v<T<Patches>> && ...),
            "The MultipatchField type should only contain instances of objects that can be "
            "manipulated like fields.");

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
    friend class MultipatchField;

    static_assert(
            !is_mem_type_v<typename base_type::example_element>,
            "For correct GPU handling a FieldMem object must be saved in a MultipatchFieldMem "
            "type.");

public:
    using span_type = MultipatchField<InternalFieldOnPatch, Patches...>;
    using view_type = MultipatchField<InternalConstFieldOnPatch, Patches...>;
    using discrete_domain_type = MultipatchType<InternalIdxRangeOnPatch, Patches...>;
    using memory_space = typename base_type::example_element::memory_space;
    using element_type = typename base_type::example_element::element_type;

    template <class Patch>
    using idx_type = typename InternalIdxRangeOnPatch<Patch>::discrete_element_type;

public:
    explicit KOKKOS_FUNCTION MultipatchField(T<Patches>... args) : base_type(args...) {}

    template <class MultipatchObj, std::enable_if_t<!is_mem_type_v<MultipatchObj>, bool> = true>
    KOKKOS_FUNCTION MultipatchField(MultipatchObj& other)
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    template <class MultipatchObj, std::enable_if_t<is_mem_type_v<MultipatchObj>, bool> = true>
    explicit MultipatchField(MultipatchObj& other)
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    template <template <typename P> typename OtherType, class... OPatches>
    MultipatchField(MultipatchField<OtherType, OPatches...>&& other)
        : base_type(std::make_tuple(other.template get<Patches>()...))
    {
        static_assert(
                std::is_same_v<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<OPatches...>>,
                "Cannot create a MultipatchField from a temporary MultipatchField with a different "
                "ordering");
        static_assert(
                std::is_same_v<std::tuple<T<Patches>...>, std::tuple<OtherType<OPatches>...>>,
                "MultipatchFields are not equivalent");
    }

    KOKKOS_DEFAULTED_FUNCTION ~MultipatchField() noexcept = default;

    template <class Patch>
    KOKKOS_FUNCTION auto get() const
    {
        return ::get_field(std::get<T<Patch>>(base_type::m_tuple));
    }

    auto idx_range() const
    {
        return MultipatchType<InternalIdxRangeOnPatch, Patches...>(
                get_idx_range(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    KOKKOS_FUNCTION auto get_field()
    {
        return MultipatchField<InternalFieldOnPatch, Patches...>(
                ::get_field(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    KOKKOS_FUNCTION auto span_view()
    {
        return get_field();
    }

    KOKKOS_FUNCTION auto get_const_field() const
    {
        return MultipatchField<InternalConstFieldOnPatch, Patches...>(
                ::get_const_field(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    KOKKOS_FUNCTION auto span_cview()
    {
        return get_const_field();
    }
};

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_multipatch_type<MultipatchField<T, Patches...>> = true;

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_data_access_methods<MultipatchField<T, Patches...>> = true;

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_multipatch_field<MultipatchField<T, Patches...>> = true;

namespace ddcHelper {

template <template <typename P> typename T1, template <typename P> typename T2, class... Patches>
void deepcopy(MultipatchField<T1, Patches...> dst, MultipatchField<T2, Patches...> src)
{
    if constexpr (ddc::is_chunk_v<typename MultipatchField<T1, Patches...>::example_element>) {
        (ddc::parallel_deepcopy(dst.template get<Patches>(), src.template get<Patches>()), ...);
    } else {
        (ddcHelper::deepcopy(dst.template get<Patches>(), src.template get<Patches>()), ...);
    }
}

} // namespace ddcHelper
```


