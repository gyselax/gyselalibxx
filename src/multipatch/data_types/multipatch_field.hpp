// SPDX-License-Identifier: MIT

#pragma once
#include "multipatch_type.hpp"


template <class T>
inline constexpr bool enable_multipatch_field = false;

template <class T>
inline constexpr bool is_multipatch_field_v
        = enable_multipatch_field<std::remove_const_t<std::remove_reference_t<T>>>;

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

    /// An internal type alias that is only instantiated if the idx_range method is called.
    template <class Patch>
    using InternalIdxRangeOnPatch = typename T<Patch>::discrete_domain_type;

    /// An internal type alias that is only instantiated if the get_const_field method is called.
    template <class Patch>
    using InternalFieldOnPatch = typename T<Patch>::span_type;

    /// An internal type alias that is only instantiated if the get_const_field method is called.
    template <class Patch>
    using InternalConstFieldOnPatch = typename T<Patch>::view_type;

    template <template <typename P> typename OtherType, class... OPatches>
    friend class MultipatchField;

    static_assert(
            !is_mem_type_v<typename base_type::example_element>,
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
    using memory_space = typename base_type::example_element::memory_space;
    /// The type of the elements inside the field.
    using element_type = typename base_type::example_element::element_type;

    /// The type used to index the field on the specified patch.
    template <class Patch>
    using idx_type = typename InternalIdxRangeOnPatch<Patch>::discrete_element_type;

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
     * This function is not explicit as it is helpful to be able to change between equivalent multipatch
     * definitions if the internal type is the same but the definition comes from different locations in
     * the code.
     *
     * @param other The equivalent MultipatchField being copied.
     */
    template <class MultipatchObj, std::enable_if_t<!is_mem_type_v<MultipatchObj>, bool> = true>
    KOKKOS_FUNCTION MultipatchField(MultipatchObj& other)
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
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
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    /**
     * Create a MultipatchField class from an r-value (temporary) instance of another MultipatchField which
     * uses the same type for the internal tuple.
     *
     * @param other The equivalent MultipatchField being copied.
     */
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

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    KOKKOS_FUNCTION auto get() const
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
     * @brief Get a MultipatchField containing modifiable fields.
     * This function matches the DDC name to allow the global get_const_field to be defined.
     *
     * @returns A set of modifiable fields providing access to the fields stored in this class.
     */
    KOKKOS_FUNCTION auto span_view()
    {
        return get_field();
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

    /**
     * @brief Get a MultipatchField containing constant fields so the values cannot be modified.
     * This function matches the DDC name to allow the global get_const_field to be defined.
     *
     * @returns A set of constant fields providing access to the fields stored in this class.
     */
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

/**
 * @brief Copy the data from one MultipatchField into another
 * @param dst The MultipatchField that the data will be copied to.
 * @param src The MultipatchField that the data will be copied from.
 */
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
