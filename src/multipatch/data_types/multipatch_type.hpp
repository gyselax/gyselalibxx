// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "types.hpp"


/**
 * @brief A class to store several objects that are of a type which is templated by the patch.
 * 
 * On a multipatch domain when we have objects and types defined on different patches, e.g. fields.
 * They can be stored in this class and then be accessed by the patch they are defined
 * on.
 * 
 * @tparam T The type of the objects that are stored on the given patches.
 * @tparam Patches The patches of the objects in the same order of the patches 
 *                 that the given objects are defined on. 
 *         
 * @warning The objects have to be defined on different patches. Otherwise retrieving 
 *          them by their patch is ill-defined.
 */
template <template <typename P> typename T, class... Patches>
class MultipatchType
{
public:
    /// @brief A tag storing the order of Patches in this MultipatchType
    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;

private:
    std::tuple<T<Patches>...> m_tuple;

    /// An internal type alias that is only instantiated if the idx_range method is called.
    template <class Patch>
    using InternalIdxRangeOnPatch = typename T<Patch>::discrete_domain_type;

    /// An internal type alias that is only instantiated if the get_const_field method is called.
    template <class Patch>
    using InternalConstFieldOnPatch = typename T<Patch>::view_type;

    template <template <typename P> typename OT, class... OPatches>
    friend class MultipatchType;

public:
    /**
     * Instantiate the MultipatchType class from an arbitrary number of objects.
     * 
     * @param args The objects to be stored in the class.
     */
    explicit MultipatchType(T<Patches>... args) : m_tuple(std::make_tuple(args...)) {}

    /**
     * Create a MultipatchType class by copying an instance of another compatible MultipatchType.
     *
     * A compatible MultipatchType is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchType. Further the original
     * MultipatchType must store objects of the correct type (the type template may be different
     * but return the same type depending on how it is designed.
     * 
     * @param other The equivalent MultipatchType being copied.
     */
    template <template <typename P> typename OT, class... OPatches>
    KOKKOS_FUNCTION MultipatchType(MultipatchType<OT, OPatches...> const& other)
        : m_tuple(std::make_tuple(other.template get<Patches>()...))
    {
        static_assert(
                ddc::type_seq_contains_v<PatchOrdering, ddc::detail::TypeSeq<OPatches...>>,
                "The type being copied does not contain all the required patches");
        static_assert(
                std::is_same_v<std::tuple<T<Patches>...>, std::tuple<OT<Patches>...>>,
                "MultipatchTypes are not equivalent");
    }

    /**
     * Create a MultipatchType class from an r-value (temporary) instance of another MultipatchType which
     * uses the same type for the internal tuple.
     * 
     * @param other The equivalent MultipatchType being copied.
     */
    template <template <typename P> typename OT, class... OPatches>
    MultipatchType(MultipatchType<OT, OPatches...>&& other) : m_tuple(std::move(other.m_tuple))
    {
        static_assert(
                std::is_same_v<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<OPatches...>>,
                "Cannot create a MultipatchType from a temporary MultipatchType with a different "
                "ordering");
        static_assert(
                std::is_same_v<std::tuple<T<Patches>...>, std::tuple<OT<OPatches>...>>,
                "MultipatchTypes are not equivalent");
    }

    KOKKOS_FUNCTION ~MultipatchType() {}

    /**
     * Retrieve an object from the patch that it is defined on.
     * 
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    KOKKOS_FUNCTION T<Patch> get() const
    {
        return std::get<T<Patch>>(m_tuple);
    }

    /**
     * @brief Get the number of objects stored inside the class. This is equal to the number of patches.
     * @return Number of elements stored in the tuple of the class.
     */
    std::size_t size() const
    {
        return sizeof...(Patches);
    }

    /**
     * @brief Get a MultipatchType containing the index ranges on which the fields are defined.
     * This function can only be called if the current MultipatchType stores fields.
     *
     * @returns The set of index ranges on which the set of fields stored in this class are defined.
     */
    auto idx_range() const
    {
        return MultipatchType<InternalIdxRangeOnPatch, Patches...>(
                get_idx_range(std::get<T<Patches>>(m_tuple))...);
    }

    /**
     * @brief Get a MultipatchType containing constant fields so the values cannot be modified.
     * This function can only be called if the current MultipatchType stores fields.
     *
     * @returns A set of constant fields providing access to the fields stored in this class.
     */
    auto get_const_field() const
    {
        return MultipatchType<InternalConstFieldOnPatch, Patches...>(
                ::get_const_field(std::get<T<Patches>>(m_tuple))...);
    }

    /**
     * @brief Get a constant reference to the tuple of objects stored inside this MultipatchType.
     *
     * @returns A constant reference to the tuple of objects stored inside this MultipatchType.
     */
    std::tuple<T<Patches>...> const& get_tuple() const
    {
        return m_tuple;
    }
};
