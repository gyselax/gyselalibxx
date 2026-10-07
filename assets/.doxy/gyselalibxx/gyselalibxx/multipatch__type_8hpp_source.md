

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

namespace detail {

template <class TypeSeqPatches, class TypeSeqInternalTypes>
class MultipatchType;

template <class... Patches, class... InternalTypes>
class MultipatchType<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<InternalTypes...>>
{
    static_assert(
            sizeof...(Patches) == sizeof...(InternalTypes),
            "There must be one internal type per patch");

public:
    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;

private:
    using InternalTypeOrdering = ddc::detail::TypeSeq<InternalTypes...>;

public:
    template <class Patch>
    using TypeOnPatch = ddc::
            type_seq_element_t<ddc::type_seq_rank_v<Patch, PatchOrdering>, InternalTypeOrdering>;

    using example_element = ddc::type_seq_element_t<0, InternalTypeOrdering>;

protected:
    std::tuple<InternalTypes...> m_tuple;

    KOKKOS_FUNCTION explicit MultipatchType(std::tuple<InternalTypes...>&& tuple) : m_tuple(tuple)
    {
    }

public:
    explicit KOKKOS_FUNCTION MultipatchType(InternalTypes... args) : m_tuple(std::move(args)...) {}

    template <class OPatchSeq, class OTypeSeq>
    KOKKOS_FUNCTION MultipatchType(MultipatchType<OPatchSeq, OTypeSeq> const& other)
        : m_tuple(other.template get<Patches>()...)
    {
        static_assert(
                ddc::type_seq_contains_v<PatchOrdering, OPatchSeq>,
                "The type being copied does not contain all the required patches");
        static_assert(
                std::is_same_v<
                        InternalTypeOrdering,
                        ddc::detail::TypeSeq<decltype(other.template get<Patches>())...>>,
                "MultipatchTypes are not equivalent");
    }

    KOKKOS_DEFAULTED_FUNCTION ~MultipatchType() noexcept = default;

    template <class Patch>
    KOKKOS_FUNCTION TypeOnPatch<Patch> get() const
            requires(!has_data_access_methods_v<TypeOnPatch<Patch>>)
    {
        return std::get<TypeOnPatch<Patch>>(m_tuple);
    }

    static constexpr std::size_t size()
    {
        return sizeof...(Patches);
    }

    KOKKOS_FUNCTION std::tuple<InternalTypes...> const& get_tuple() const
    {
        return m_tuple;
    }
};

} // namespace detail

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool
        enable_multipatch_type<detail::MultipatchType<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <template <typename P> typename T, class... Patches>
using MultipatchType = detail::
        MultipatchType<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...>>;
```


