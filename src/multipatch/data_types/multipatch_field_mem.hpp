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

/**
 * @brief A class to store field memory block objects on patches.
 *
 * On a multipatch domain when we have objects and types defined on different patches, e.g. FieldMems.
 * They can be stored in this class and then be accessed by the patch they are defined
 * on.
 *
 * This class should not be used directly to instantiate objects. Instead the global alias
 * MultipatchFieldMem should be used.
 * I.e. use `MultipatchFieldMem<TypeOnPatch, Patches...>` instead of
 * `detail::MultipatchFieldMem<TypeSeqPatches, TypeSeqInternalTypes>`.
 *
 * @tparam Patches The patches on which the objects are defined.
 * @tparam InternalTypes The types of the memory allocating field objects that are stored on the given patches.
 *                 The order must match that of the patches.
 *
 * @warning The objects have to be defined on different patches. Otherwise retrieving
 *          them by their patch is ill-defined.
 */
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
    /// @brief The MultipatchType from which this class inherits
    using base_type = MultipatchType<
            ddc::detail::TypeSeq<Patches...>,
            ddc::detail::TypeSeq<InternalTypes...>>;

    /// @brief A tag storing the order of Patches in this MultipatchFieldMem
    using typename base_type::PatchOrdering;

    /// @brief The type of the object stored on the given patch.
    template <class Patch>
    using TypeOnPatch = typename base_type::template TypeOnPatch<Patch>;

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
     * Instantiate the MultipatchFieldMem class from an arbitrary number of objects.
     *
     * @param label A label used to tag parallel regions and memory allocations for profiling.
     * @param args The objects to be stored in the class.
     */
    explicit MultipatchFieldMem(std::string const& label, InternalTypes... args)
        : base_type(label, args...)
    {
    }

    /// Version without a label
    explicit MultipatchFieldMem(InternalTypes... args) : MultipatchFieldMem("no-label", args...) {}

    /**
     * Create a MultipatchFieldMem class by copying an instance of another compatible MultipatchFieldMem.
     *
     * A compatible MultipatchFieldMem is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchFieldMem.
     *
     * @param label A label used to tag parallel regions and memory allocations for profiling.
     * @param other The equivalent MultipatchFieldMem being copied.
     */
    template <class MultipatchObj>
    explicit MultipatchFieldMem(std::string const& label, MultipatchObj& other)
        : base_type(InternalTypes(label, other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    /// Version without a label
    template <class MultipatchObj>
    explicit MultipatchFieldMem(MultipatchObj& other) : MultipatchFieldMem("no-label", other)
    {
    }

    ~MultipatchFieldMem() noexcept = default;

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    typename TypeOnPatch<Patch>::view_type get() const
    {
        return ::get_const_field(std::get<TypeOnPatch<Patch>>(base_type::m_tuple));
    }

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    typename TypeOnPatch<Patch>::span_type get()
    {
        return ::get_field(std::get<TypeOnPatch<Patch>>(base_type::m_tuple));
    }

    /**
     * @brief Get a MultipatchType containing the index ranges on which the fields are defined.
     *
     * @returns The set of index ranges on which the set of fields stored in this class are defined.
     */
    discrete_domain_type idx_range() const
    {
        return discrete_domain_type(
                get_idx_range(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    /**
     * @brief Get a MultipatchField containing modifiable fields.
     *
     * @returns A set of modifiable fields providing access to the fields stored in this class.
     */
    span_type get_field()
    {
        return span_type(::get_field(std::get<TypeOnPatch<Patches>>(base_type::m_tuple))...);
    }

    /**
     * @brief Get a MultipatchConstField containing constant fields so the values cannot be modified.
     *
     * @returns A set of constant fields providing access to the fields stored in this class.
     */
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

/**
 * @brief A class to store field memory block objects on patches.
 *
 * See detail::MultipatchFieldMem for more details. This alias creates the type T<Patch> for each
 * of the patches. Two type templates which return the same type for each patch therefore lead to
 * the same MultipatchFieldMem.
 *
 * @tparam T The type of the FieldMem/DerivMem/VectorFieldMem that are stored on the given patches.
 * @tparam Patches The patches of the objects in the same order of the patches
 *                 that the given objects are defined on.
 */
template <template <typename P> typename T, class... Patches>
using MultipatchFieldMem = detail::
        MultipatchFieldMem<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...>>;
