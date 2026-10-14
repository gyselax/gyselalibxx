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

/**
 * @brief A class to store several objects that are of a type which is templated by the patch.
 *
 * On a multipatch domain when we have objects and types defined on different patches, e.g. fields.
 * They can be stored in this class and then be accessed by the patch they are defined
 * on.
 *
 * This class should not be used directly for instantiation. Instead the global alias MultipatchType
 * should be used. This alias packs the patches and the types into TypeSeqs. This ensures that
 * equivalent type templates (which return the same type for each patch) lead to the same class.
 *
 * I.e. use `MultipatchType<TypeOnPatch, Patches...>` instead of
 * `detail::MultipatchType<TypeSeqPatches, TypeSeqInternalTypes>`.
 *
 * @tparam Patches The patches on which the objects are defined.
 * @tparam InternalTypes The types of the objects that are stored on the given patches.
 *                 The order must match that of the patches.
 *
 * @warning The objects have to be defined on different patches. Otherwise retrieving
 *          them by their patch is ill-defined.
 */
template <class... Patches, class... InternalTypes>
class MultipatchType<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<InternalTypes...>>
{
    static_assert(
            sizeof...(Patches) == sizeof...(InternalTypes),
            "There must be one internal type per patch");

public:
    /// @brief A tag storing the order of Patches in this MultipatchType
    using PatchOrdering = ddc::detail::TypeSeq<Patches...>;

private:
    /// @brief A tag storing the order of the internal types in this MultipatchType
    using InternalTypeOrdering = ddc::detail::TypeSeq<InternalTypes...>;

public:
    /// @brief The type of the object stored on the given patch.
    template <class Patch>
    using TypeOnPatch = ddc::
            type_seq_element_t<ddc::type_seq_rank_v<Patch, PatchOrdering>, InternalTypeOrdering>;

    /**
     * @brief The type of one of the elements of the MultipatchType. This can be used to check that
     * types are as expected using functions such as ddc::is_chunk_v.
     */
    using example_element = ddc::type_seq_element_t<0, InternalTypeOrdering>;

protected:
    /// The internal tuple containing the data
    std::tuple<InternalTypes...> m_tuple;

    /**
     * A constructor for sub-classes which can build the necessary tuple directly following their own rules.
     *
     * @param tuple The internal tuple.
     */
    KOKKOS_FUNCTION explicit MultipatchType(std::tuple<InternalTypes...>&& tuple) : m_tuple(tuple)
    {
    }

public:
    /**
     * Instantiate the MultipatchType class from an arbitrary number of objects.
     *
     * @param args The objects to be stored in the class.
     */
    explicit KOKKOS_FUNCTION MultipatchType(InternalTypes... args) : m_tuple(std::move(args)...) {}

    /**
     * Create a MultipatchType class by copying an instance of another compatible MultipatchType.
     *
     * A compatible MultipatchType is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchType. Further the original
     * MultipatchType must store objects of the correct type (the type template may be different
     * but return the same type depending on how it is designed.
     *
     * This function is not explicit as it is helpful to be able to change between equivalent multipatch
     * definitions if the internal type is the same but the definition comes from different locations in
     * the code.
     *
     * @param other The equivalent MultipatchType being copied.
     */
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

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    KOKKOS_FUNCTION TypeOnPatch<Patch> get() const
            requires(!has_data_access_methods_v<TypeOnPatch<Patch>>)
    {
        return std::get<TypeOnPatch<Patch>>(m_tuple);
    }

    /**
     * @brief Get the number of objects stored inside the class. This is equal to the number of patches.
     * @return Number of elements stored in the tuple of the class.
     */
    static constexpr std::size_t size()
    {
        return sizeof...(Patches);
    }

    /**
     * @brief Get a constant reference to the tuple of objects stored inside this MultipatchType.
     *
     * @returns A constant reference to the tuple of objects stored inside this MultipatchType.
     */
    KOKKOS_FUNCTION std::tuple<InternalTypes...> const& get_tuple() const
    {
        return m_tuple;
    }
};

} // namespace detail

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool
        enable_multipatch_type<detail::MultipatchType<TypeSeqPatches, TypeSeqInternalTypes>> = true;

/**
 * @brief A class to store several objects that are of a type which is templated by the patch.
 *
 * See detail::MultipatchType for more details. This alias creates the type T<Patch> for each of
 * the patches. Two type templates which return the same type for each patch therefore lead to the
 * same MultipatchType.
 *
 * @tparam T The type of the objects that are stored on the given patches.
 * @tparam Patches The patches of the objects in the same order of the patches
 *                 that the given objects are defined on.
 */
template <template <typename P> typename T, class... Patches>
using MultipatchType = detail::
        MultipatchType<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...>>;
