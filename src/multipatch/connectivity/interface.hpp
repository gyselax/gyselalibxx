// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "edge.hpp"

/**
 * @brief Represent a simple sticking of two edges.
 * 
 * @tparam Edge1Type Edge type defined on first patch.
 * @tparam Edge2Type Edge type defined on second patch.
 * @tparam orientations_agree_bool Boolean, 
 *      if true, the parametrisations of the physical edge have the same orientation; 
 *      else, the parametrisations of the physical edge have the opposite orientation.
 * @see EdgeTransformation
*/
template <class Edge1Type, class Edge2Type, bool orientations_agree_bool>
struct Interface
{
    static_assert(
            !std::is_same_v<Edge1Type, OutsideEdge> || !std::is_same_v<Edge2Type, OutsideEdge>,
            "At least one edge should be inside the simulation.");

    /// @brief Edge type of the first patch.
    using Edge1 = Edge1Type;
    /// @brief Edge type of the second patch.
    using Edge2 = Edge2Type;

    /**
     * If True, the parametrisations of the physical edge have the same orientation.
     * If False, the parametrisations of the physical edge have the opposite orientation.
     * (See EdgeTransformation).
    */
    static constexpr bool orientations_agree = orientations_agree_bool;

    template <class Patch>
    static constexpr bool connected_to_patch()
    {
        if constexpr (std::is_same_v<Edge1, OutsideEdge>) {
            return std::is_same_v<typename Edge2::associated_patch, Patch>;
        } else if constexpr (std::is_same_v<Edge2, OutsideEdge>) {
            return std::is_same_v<typename Edge1::associated_patch, Patch>;
        } else {
            return (std::is_same_v<typename Edge1::associated_patch, Patch>)
                   || (std::is_same_v<typename Edge2::associated_patch, Patch>);
        }
    }
};
