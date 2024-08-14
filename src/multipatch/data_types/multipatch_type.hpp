// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"


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
    std::tuple<T<Patches>...> m_tuple;

public:
    /**
     * Instantiate the MultipatchType class from an arbitrary number of objects.
     * 
     * @param args The objects to be stored in the class.
     */
    MultipatchType(T<Patches>... args) : m_tuple(std::make_tuple(args...)) {}

    /**
     * Retrieve an object from the patch that it is defined on.
     * 
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class Patch>
    T<Patch> get() const
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
};
