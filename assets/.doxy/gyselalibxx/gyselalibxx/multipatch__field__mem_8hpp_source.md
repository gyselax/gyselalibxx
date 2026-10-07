

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

namespace detail {

template <class TypeSeqPatches, class TypeSeqInternalTypes>
class MultipatchFieldMem;

template <class... Patches, class... InternalTypes>
class MultipatchFieldMem<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<InternalTypes...>>
    : public MultipatchType<
              ddc::detail::TypeSeq<Patches...>,
              ddc::detail::TypeSeq<InternalTypes...>>
{
    static_assert(
            (has_data_access_methods_v<InternalTypes> && ...),
            "The MultipatchFieldMem type should only contain instances of objects that can be "
            "manipulated like fields.");
    static_assert(
            (is_mem_type_v<InternalTypes> && ...),
            "The MultipatchFieldMem type should only contain instances of objects that allocate "
            "memory.");

public:
    using base_type = MultipatchType<
            ddc::detail::TypeSeq<Patches...>,
            ddc::detail::TypeSeq<InternalTypes...>>;

    using typename base_type::PatchOrdering;

    template <class Patch>
    using TypeOnPatch = typename base_type::template TypeOnPatch<Patch>;

public:
    using span_type = MultipatchField<
            PatchOrdering,
            ddc::detail::TypeSeq<typename InternalTypes::span_type...>>;
    using view_type = MultipatchField<
            PatchOrdering,
            ddc::detail::TypeSeq<typename InternalTypes::view_type...>>;
    using discrete_domain_type = MultipatchType<
            PatchOrdering,
            ddc::detail::TypeSeq<typename InternalTypes::discrete_domain_type...>>;
    using memory_space = typename base_type::example_element::memory_space;
    using element_type = typename base_type::example_element::element_type;

public:
    explicit MultipatchFieldMem(std::string const& label, InternalTypes... args)
        : base_type(label, args...)
    {
    }

    explicit MultipatchFieldMem(InternalTypes... args) : MultipatchFieldMem("no-label", args...) {}

    template <class MultipatchObj>
    explicit MultipatchFieldMem(std::string const& label, MultipatchObj& other)
        : base_type(InternalTypes(label, other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    template <class MultipatchObj>
    explicit MultipatchFieldMem(MultipatchObj& other) : MultipatchFieldMem("no-label", other)
    {
    }

    ~MultipatchFieldMem() noexcept = default;

    template <class Patch>
    typename TypeOnPatch<Patch>::view_type get() const
    {
        return ::get_const_field(std::get<TypeOnPatch<Patch>>(base_type::m_tuple));
    }

    template <class Patch>
    typename TypeOnPatch<Patch>::span_type get()
    {
        return ::get_field(std::get<TypeOnPatch<Patch>>(base_type::m_tuple));
    }

    discrete_domain_type idx_range() const
    {
        return discrete_domain_type(
                get_idx_range(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    span_type get_field()
    {
        return span_type(::get_field(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    view_type get_const_field() const
    {
        return view_type(::get_const_field(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }
};

} // namespace detail

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_multipatch_type<
        detail::MultipatchFieldMem<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_multipatch_field_mem<
        detail::MultipatchFieldMem<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool
        enable_mem_type<detail::MultipatchFieldMem<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_data_access_methods<
        detail::MultipatchFieldMem<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <template <typename P> typename T, class... Patches>
using MultipatchFieldMem = detail::
        MultipatchFieldMem<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...>>;
```


