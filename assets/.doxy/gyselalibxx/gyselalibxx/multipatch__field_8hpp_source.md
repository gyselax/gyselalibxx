

# File multipatch\_field.hpp

[**File List**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_field.hpp**](multipatch__field_8hpp.md)

[Go to the documentation of this file](multipatch__field_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include "multipatch_type.hpp"
#include "patch.hpp"


template <class T>
inline constexpr bool enable_multipatch_field = false;

template <class T>
inline constexpr bool is_multipatch_field_v
        = enable_multipatch_field<std::remove_const_t<std::remove_reference_t<T>>>;

namespace detail {

template <class TypeSeqPatches, class TypeSeqInternalTypes>
class MultipatchField;

template <class... Patches, class... InternalTypes>
class MultipatchField<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<InternalTypes...>>
    : public MultipatchType<
              ddc::detail::TypeSeq<Patches...>,
              ddc::detail::TypeSeq<InternalTypes...>>
{
    static_assert(
            (has_data_access_methods_v<InternalTypes> && ...),
            "The MultipatchField type should only contain instances of objects that can be "
            "manipulated like fields.");

public:
    using base_type = MultipatchType<
            ddc::detail::TypeSeq<Patches...>,
            ddc::detail::TypeSeq<InternalTypes...>>;

    using typename base_type::PatchOrdering;

    template <class Patch>
    using TypeOnPatch = typename base_type::template TypeOnPatch<Patch>;

    static_assert(
            !is_mem_type_v<typename base_type::example_element>,
            "For correct GPU handling a FieldMem object must be saved in a MultipatchFieldMem "
            "type.");

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
    explicit KOKKOS_FUNCTION MultipatchField(InternalTypes... args) : base_type(args...) {}

    template <class MultipatchObj, std::enable_if_t<!is_mem_type_v<MultipatchObj>, bool> = true>
    explicit KOKKOS_FUNCTION MultipatchField(MultipatchObj& other)
        : base_type(InternalTypes(other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    template <class MultipatchObj, std::enable_if_t<is_mem_type_v<MultipatchObj>, bool> = true>
    explicit MultipatchField(MultipatchObj& other)
        : base_type(InternalTypes(other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    KOKKOS_DEFAULTED_FUNCTION MultipatchField(MultipatchField const&) noexcept = default;

    KOKKOS_DEFAULTED_FUNCTION MultipatchField(MultipatchField&&) noexcept = default;

    KOKKOS_DEFAULTED_FUNCTION ~MultipatchField() noexcept = default;

    template <class PatchType>
    KOKKOS_FUNCTION auto get() const requires(is_patch_v<PatchType>)
    {
        return std::get<TypeOnPatch<PatchType>>(base_type::m_tuple);
    }

    KOKKOS_FUNCTION discrete_domain_type idx_range() const
    {
        return discrete_domain_type(
                get_idx_range(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    KOKKOS_FUNCTION span_type get_field()
    {
        return span_type(::get_field(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    KOKKOS_FUNCTION view_type get_const_field() const
    {
        return view_type(::get_const_field(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    template <class QueryTag>
    inline constexpr auto get() const noexcept requires(
            (is_vector_field_v<InternalTypes> && ...)
            && (ddc::in_tags_v<QueryTag, typename InternalTypes::NDTypeTag> && ...))
    {
        using FieldType = MultipatchField<
                PatchOrdering,
                ddc::detail::TypeSeq<typename InternalTypes::chunk_span_type...>>;
        return FieldType(ddcHelper::get<QueryTag>(std::get<InternalTypes>(base_type::m_tuple))...);
    }
};

} // namespace detail

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_multipatch_type<
        detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_data_access_methods<
        detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_multipatch_field<
        detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <template <typename P> typename T, class... Patches>
using MultipatchField = detail::
        MultipatchField<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...>>;

template <class... VectorFieldType, class PatchTypeSeq>
inline constexpr bool enable_vector_field<detail::MultipatchField<
        PatchTypeSeq,
        ddc::detail::TypeSeq<VectorFieldType...>>> = (is_vector_field_v<VectorFieldType> && ...);

namespace ddcHelper {

template <class... Patches, class TypeSeqInternalTypes1, class TypeSeqInternalTypes2>
void deepcopy(
        detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes1> dst,
        detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes2> src)
{
    using DstType
            = detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes1>;
    if constexpr (ddc::is_chunk_v<typename DstType::example_element>) {
        (ddc::parallel_deepcopy(dst.template get<Patches>(), src.template get<Patches>()), ...);
    } else {
        (ddcHelper::deepcopy(dst.template get<Patches>(), src.template get<Patches>()), ...);
    }
}

} // namespace ddcHelper
```


