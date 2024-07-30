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

    /// @brief Type of first patch.
    using Patch1 = typename Edge1::associated_patch;
    /// @brief Type of second patch.
    using Patch2 = typename Edge2::associated_patch;
};
