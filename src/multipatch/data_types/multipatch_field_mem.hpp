// SPDX-License-Identifier: MIT
#pragma once

#include "multipatch_field.hpp"
#include "multipatch_type.hpp"

template <class T>
inline constexpr bool enable_multipatch_field_mem = false;

template <class T>
inline constexpr bool is_multipatch_field_mem_v
        = enable_multipatch_field_mem<std::remove_const_t<std::remove_reference_t<T>>>;

/**
 * @brief A class to store field memory block objects on patches.
 *
 * On a multipatch domain when we have objects and types defined on different patches, e.g. FieldMems.
 * They can be stored in this class and then be accessed by the patch they are defined
 * on.
 *
 * @tparam T The type of the FieldMem/DerivMem/VectorFieldMem that are stored on the given patches.
 * @tparam Patches The patches of the objects in the same order of the patches
 *                 that the given objects are defined on.
 *
 * @warning The objects have to be defined on different patches. Otherwise retrieving
 *          them by their patch is ill-defined.
 */
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
    /// @brief The MultipatchType from which this class inherits
    using base_type = MultipatchType<T, Patches...>;

    /// @brief A tag storing the order of Patches in this MultipatchFieldMem
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
    friend class MultipatchFieldMem;

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

public:
    /**
     * Instantiate the MultipatchFieldMem class from an arbitrary number of objects.
     *
     * @param args The objects to be stored in the class.
     */
    explicit MultipatchFieldMem(T<Patches>... args) : base_type(args...) {}

    /**
     * Create a MultipatchFieldMem class by copying an instance of another compatible MultipatchFieldMem.
     *
     * A compatible MultipatchFieldMem is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchFieldMem. Further the original
     * MultipatchFieldMem must store objects of the correct type (the type template may be different
     * but return the same type depending on how it is designed.
     *
     * @param other The equivalent MultipatchFieldMem being copied.
     */
    template <class MultipatchObj>
    explicit MultipatchFieldMem(MultipatchObj& other)
        : base_type(std::move(T<Patches>(other.template get<Patches>()))...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    /**
     * Create a MultipatchFieldMem class from an r-value (temporary) instance of another MultipatchFieldMem which
     * uses the same type for the internal tuple.
     *
     * @param other The equivalent MultipatchFieldMem being copied.
     */
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

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    auto get() const
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
    auto get()
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
    auto get_field()
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
    auto span_view()
    {
        return get_field();
    }

    /**
     * @brief Get a MultipatchConstField containing constant fields so the values cannot be modified.
     *
     * @returns A set of constant fields providing access to the fields stored in this class.
     */
    auto get_const_field() const
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
