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

protected:
    /// The internal tuple containing the data
    std::tuple<T<Patches>...> m_tuple;

    template <template <typename P> typename OT, class... OPatches>
    friend class MultipatchType;

public:
    /**
     * Instantiate the MultipatchType class from an arbitrary number of objects.
     * 
     * @param args The objects to be stored in the class.
     */
    explicit KOKKOS_FUNCTION MultipatchType(T<Patches>... args)
        : m_tuple(std::make_tuple(std::move(args)...))
    {
    }

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
     * @brief Get a constant reference to the tuple of objects stored inside this MultipatchType.
     *
     * @returns A constant reference to the tuple of objects stored inside this MultipatchType.
     */
    KOKKOS_FUNCTION std::tuple<T<Patches>...> const& get_tuple() const
    {
        return m_tuple;
    }
};

template <template <typename P> typename T, class... Patches>
inline constexpr bool enable_multipatch_type<MultipatchType<T, Patches...>> = true;
