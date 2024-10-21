// SPDX-License-Identifier: MIT

#pragma once
#include "multipatch_type.hpp"


/**
 * @brief A class to store field objects on patches.
 *
 * On a multipatch domain when we have objects and types defined on different patches, e.g. fields.
 * They can be stored in this class and then be accessed by the patch they are defined
 * on.
 *
 * @tparam T The type of the fields/derivative fields/vector fields that are stored on the given patches.
 * @tparam Patches The patches of the objects in the same order of the patches
 *                 that the given objects are defined on.
 *
 * @warning The objects have to be defined on different patches. Otherwise retrieving
 *          them by their patch is ill-defined.
 */
template <template <typename P> typename T, class... Patches>
class MultipatchField : public MultipatchType<T, Patches...>
{
    static_assert(
            (has_data_access_methods_v<T<Patches>> && ...),
            "The MultipatchField type should only contain instances of objects that can be "
            "manipulated like fields.");

public:
    /// @brief The MultipatchType from which this class inherits
    using base_type = MultipatchType<T, Patches...>;

    /// @brief A tag storing the order of Patches in this MultipatchField
    using typename base_type::PatchOrdering;

private:
    /// An internal type alias that is only instantiated if the idx_range method is called.
    template <class Patch>
    using InternalIdxRangeOnPatch = typename T<Patch>::discrete_domain_type;

    /// An internal type alias that is only instantiated if the get_const_field method is called.
    template <class Patch>
    using InternalFieldOnPatch = typename T<Patch>::span_type;

    /// An internal type alias that is only instantiated if the get_const_field method is called.
    template <class Patch>
    using InternalConstFieldOnPatch = typename T<Patch>::view_type;

    template <template <typename P> typename OT, class... OPatches>
    friend class MultipatchField;

    using example_element = T<ddc::type_seq_element_t<0, PatchOrdering>>;

    static_assert(
            !is_mem_type_v<example_element>,
            "For correct GPU handling a FieldMem object must be saved in a MultipatchFieldMem "
            "type.");

public:
    /// The type of a modifiable reference to this multipatch field
    using span_type = MultipatchField<InternalFieldOnPatch, Patches...>;
    /// The type of a constant reference to this multipatch field
    using view_type = MultipatchField<InternalConstFieldOnPatch, Patches...>;
    /// The type of the index ranges that can be used to access this field.
    using discrete_domain_type = MultipatchType<InternalIdxRangeOnPatch, Patches...>;
    /// The memory space (CPU/GPU) where the data is saved.
    using memory_space = typename example_element::memory_space;
    /// The type of the elements inside the field.
    using element_type = typename example_element::element_type;

public:
    /**
     * Instantiate the MultipatchField class from an arbitrary number of objects.
     *
     * @param args The objects to be stored in the class.
     */
    explicit KOKKOS_FUNCTION MultipatchField(T<Patches>... args) : base_type(args...) {}

    /**
     * Create a MultipatchField class by copying an instance of another compatible MultipatchField.
     *
     * A compatible MultipatchField is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchField. Further the original
     * MultipatchField must store objects of the correct type (the type template may be different
     * but return the same type depending on how it is designed).
     *
     * @param other The equivalent MultipatchField being copied.
     */
    template <class MultipatchObj, std::enable_if_t<!is_mem_type_v<MultipatchObj>, bool> = true>
    KOKKOS_FUNCTION MultipatchField(MultipatchObj& other)
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
    {
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
    MultipatchField(MultipatchObj& other)
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
    {
    }

    /**
     * Create a MultipatchField class from an r-value (temporary) instance of another MultipatchField which
     * uses the same type for the internal tuple.
     *
     * @param other The equivalent MultipatchField being copied.
     */
    template <template <typename P> typename OT, class... OPatches>
    MultipatchField(MultipatchField<OT, OPatches...>&& other) : base_type(other)
    {
        static_assert(
                std::is_same_v<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<OPatches...>>,
                "Cannot create a MultipatchField from a temporary MultipatchField with a different "
                "ordering");
        static_assert(
                std::is_same_v<std::tuple<T<Patches>...>, std::tuple<OT<OPatches>...>>,
                "MultipatchFields are not equivalent");
    }

    KOKKOS_FUNCTION ~MultipatchField() {}

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    KOKKOS_FUNCTION auto get() const
    {
        return ::get_const_field(std::get<T<Patch>>(base_type::m_tuple));
    }

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    KOKKOS_FUNCTION auto get()
    {
        return ::get_field(std::get<T<Patch>>(base_type::m_tuple));
    }

    /**
     * @brief Get a MultipatchType containing the index ranges on which the fields are defined.
     *
     * @returns The set of index ranges on which the set of fields stored in this class are defined.
     */
    auto idx_range() const
    {
        return MultipatchType<InternalIdxRangeOnPatch, Patches...>(
                get_idx_range(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    /**
     * @brief Get a MultipatchField containing modifiable fields.
     *
     * @returns A set of modifiable fields providing access to the fields stored in this class.
     */
    KOKKOS_FUNCTION auto get_field()
    {
        return MultipatchField<InternalFieldOnPatch, Patches...>(
                ::get_field(std::get<T<Patches>>(base_type::m_tuple))...);
    }

    /**
     * @brief Get a MultipatchField containing constant fields so the values cannot be modified.
     *
     * @returns A set of constant fields providing access to the fields stored in this class.
     */
    KOKKOS_FUNCTION auto get_const_field() const
    {
        return MultipatchField<InternalConstFieldOnPatch, Patches...>(
                ::get_const_field(std::get<T<Patches>>(base_type::m_tuple))...);
    }
};

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_multipatch_type<MultipatchField<T, Patches...>> = true;
