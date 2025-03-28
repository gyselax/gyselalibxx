

# File multipatch\_type.hpp

[**File List**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_type.hpp**](multipatch__type_8hpp.md)

[Go to the documentation of this file](multipatch__type_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "types.hpp"

template <class T>
inline constexpr bool enable_multipatch_type = false;

template <class T>
inline constexpr bool is_multipatch_type_v
        = enable_multipatch_type<std::remove_const_t<std::remove_reference_t<T>>>;


template <template <typename P> typename T, class... Patches>
class MultipatchType
{
public:
    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;

    using example_element = T<ddc::type_seq_element_t<0, PatchOrdering>>;

protected:
    std::tuple<T<Patches>...> m_tuple;

    template <template <typename P> typename OtherType, class... OPatches>
    friend class MultipatchType;

    KOKKOS_FUNCTION explicit MultipatchType(std::tuple<T<Patches>...>&& tuple) : m_tuple(tuple) {}

public:
    explicit KOKKOS_FUNCTION MultipatchType(T<Patches>... args)
        : m_tuple(std::make_tuple(std::move(args)...))
    {
    }

    template <template <typename P> typename OtherType, class... OPatches>
    KOKKOS_FUNCTION MultipatchType(MultipatchType<OtherType, OPatches...> const& other)
        : m_tuple(std::make_tuple(other.template get<Patches>()...))
    {
        static_assert(
                ddc::type_seq_contains_v<PatchOrdering, ddc::detail::TypeSeq<OPatches...>>,
                "The type being copied does not contain all the required patches");
        static_assert(
                std::is_same_v<std::tuple<T<Patches>...>, std::tuple<OtherType<Patches>...>>,
                "MultipatchTypes are not equivalent");
    }

    template <template <typename P> typename OtherType, class... OPatches>
    MultipatchType(MultipatchType<OtherType, OPatches...>&& other)
        : m_tuple(std::make_tuple(std::move(other.template get<Patches>())...))
    {
        static_assert(
                std::is_same_v<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<OPatches...>>,
                "Cannot create a MultipatchType from a temporary MultipatchType with a different "
                "ordering");
        static_assert(
                std::is_same_v<std::tuple<T<Patches>...>, std::tuple<OtherType<OPatches>...>>,
                "MultipatchTypes are not equivalent");
    }

    KOKKOS_DEFAULTED_FUNCTION ~MultipatchType() noexcept = default;

    template <class Patch, std::enable_if_t<!has_data_access_methods_v<T<Patch>>, bool> = true>
    KOKKOS_FUNCTION T<Patch> get() const
    {
        return std::get<T<Patch>>(m_tuple);
    }

    static constexpr std::size_t size()
    {
        return sizeof...(Patches);
    }

    KOKKOS_FUNCTION std::tuple<T<Patches>...> const& get_tuple() const
    {
        return m_tuple;
    }
};

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_multipatch_type<MultipatchType<T, Patches...>> = true;
```


