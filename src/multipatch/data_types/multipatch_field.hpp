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

/**
 * @brief A class to store field objects on patches.
 *
 * On a multipatch domain when we have objects and types defined on different patches, e.g. fields.
 * They can be stored in this class and then be accessed by the patch they are defined
 * on.
 *
 * This class should not be used directly to instantiate objects. Instead the global alias
 * MultipatchField should be used.
 * I.e. use `MultipatchField<TypeOnPatch, Patches...>` instead of
 * `detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>`.
 *
 * @tparam Patches The patches on which the objects are defined.
 * @tparam InternalTypes The types of the field objects that are stored on the given patches.
 *                 The order must match that of the patches.
 *
 * @warning The objects have to be defined on different patches. Otherwise retrieving
 *          them by their patch is ill-defined.
 */
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
    /// @brief The MultipatchType from which this class inherits
    using base_type = MultipatchType<
            ddc::detail::TypeSeq<Patches...>,
            ddc::detail::TypeSeq<InternalTypes...>>;

    /// @brief A tag storing the order of Patches in this MultipatchField
    using typename base_type::PatchOrdering;

    /// @brief The type of the object stored on the given patch.
    template <class Patch>
    using TypeOnPatch = typename base_type::template TypeOnPatch<Patch>;

    static_assert(
            !is_mem_type_v<typename base_type::example_element>,
            "For correct GPU handling a FieldMem object must be saved in a MultipatchFieldMem "
            "type.");

public:
    /// The type of a modifiable reference to this multipatch field
    using span_type = MultipatchField<
            PatchOrdering,
            ddc::detail::TypeSeq<typename InternalTypes::span_type...>>;
    /// The type of a constant reference to this multipatch field
    using view_type = MultipatchField<
            PatchOrdering,
            ddc::detail::TypeSeq<typename InternalTypes::view_type...>>;
    /// The type of the index ranges that can be used to access this field.
    using discrete_domain_type = MultipatchType<
            PatchOrdering,
            ddc::detail::TypeSeq<typename InternalTypes::discrete_domain_type...>>;
    /// The memory space (CPU/GPU) where the data is saved.
    using memory_space = typename base_type::example_element::memory_space;
    /// The type of the elements inside the field.
    using element_type = typename base_type::example_element::element_type;

public:
    /**
     * Instantiate the MultipatchField class from an arbitrary number of objects.
     *
     * @param args The objects to be stored in the class.
     */
    explicit KOKKOS_FUNCTION MultipatchField(InternalTypes... args) : base_type(args...) {}

    /**
     * Create a MultipatchField class by copying an instance of another compatible MultipatchField.
     *
     * A compatible MultipatchField is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchField. Further the original
     * MultipatchField must store objects which can be converted to the correct type.
     *
     * @param other The equivalent MultipatchField being copied.
     */
    template <class MultipatchObj, std::enable_if_t<!is_mem_type_v<MultipatchObj>, bool> = true>
    explicit KOKKOS_FUNCTION MultipatchField(MultipatchObj& other)
        : base_type(InternalTypes(other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    /**
     * Create a MultipatchField class from a compatible MultipatchFieldMem.
     *
     * A compatible MultipatchField is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchField. Further the original
     * MultipatchField must store objects of the correct type.
     *
     * @param other The MultipatchFieldMem being accessed.
     */
    template <class MultipatchObj, std::enable_if_t<is_mem_type_v<MultipatchObj>, bool> = true>
    explicit MultipatchField(MultipatchObj& other)
        : base_type(InternalTypes(other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    KOKKOS_DEFAULTED_FUNCTION MultipatchField(MultipatchField const&) noexcept = default;

    KOKKOS_DEFAULTED_FUNCTION MultipatchField(MultipatchField&&) noexcept = default;

    KOKKOS_DEFAULTED_FUNCTION ~MultipatchField() noexcept = default;

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class PatchType>
    KOKKOS_FUNCTION TypeOnPatch<PatchType> get() const requires(is_patch_v<PatchType>)
    {
        return ::get_field(std::get<TypeOnPatch<PatchType>>(base_type::m_tuple));
    }

    /**
     * @brief Get a MultipatchType containing the index ranges on which the fields are defined.
     *
     * @returns The set of index ranges on which the set of fields stored in this class are defined.
     */
    KOKKOS_FUNCTION discrete_domain_type idx_range() const
    {
        return discrete_domain_type(
                get_idx_range(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    /**
     * @brief Get a MultipatchField containing modifiable fields.
     *
     * @returns A set of modifiable fields providing access to the fields stored in this class.
     */
    KOKKOS_FUNCTION span_type get_field()
    {
        return span_type(::get_field(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    /**
     * @brief Get a MultipatchField containing constant fields so the values cannot be modified.
     *
     * @returns A set of constant fields providing access to the fields stored in this class.
     */
    KOKKOS_FUNCTION view_type get_const_field() const
    {
        return view_type(::get_const_field(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    /**
     * @brief Get the Field describing the component in the QueryTag direction.
     *
     * @return The field in the specified direction.
     */
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

/**
 * @brief A class to store field objects on patches.
 *
 * See detail::MultipatchField for more details. This alias creates the type T<Patch> for each of
 * the patches. Two type templates which return the same type for each patch therefore lead to the
 * same MultipatchField.
 *
 * @tparam T The type of the fields/derivative fields/vector fields that are stored on the given patches.
 * @tparam Patches The patches of the objects in the same order of the patches
 *                 that the given objects are defined on.
 */
template <template <typename P> typename T, class... Patches>
using MultipatchField = detail::
        MultipatchField<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...>>;

template <class... VectorFieldType, class PatchTypeSeq>
inline constexpr bool enable_vector_field<detail::MultipatchField<
        PatchTypeSeq,
        ddc::detail::TypeSeq<VectorFieldType...>>> = (is_vector_field_v<VectorFieldType> && ...);

namespace ddcHelper {

/**
 * @brief Copy the data from one MultipatchField into another
 * @param dst The MultipatchField that the data will be copied to.
 * @param src The MultipatchField that the data will be copied from.
 */
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
