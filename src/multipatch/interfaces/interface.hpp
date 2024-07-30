// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "edge.hpp"

/**
 * @brief Represent a simple sticking of two edges.
 * 
 * @tparam IDim1 
 *      The discrete dimension of the first edge. 
 * @tparam IDim2
 *      The discrete dimension of the second edge. 
 *  
 * @see EdgeCoordinatesTransformation
*/
template <class IDim1, class IDim2>
struct Interface
{
    /**
     * The Edge of the first patch. 
    */
    Edge<IDim1> edge_1;
    /**
     * The Edge of the second patch. 
    */
    Edge<IDim2> edge_2;
    /**
     * If True, the parametrisations of the physical edge have the same orientation.
     * If False, the parametrisations of the physical edge have the opposite orientation.
     * (See EdgeCoordinatesTransformation).
    */
    bool orientations_agree;

    /**
     * @brief Get the edge on the given dimension. 
     * 
     * @tparam IDim
     *      The discrete dimension of the edge. 
     * 
     * @return The edge of the interface defined on the given dimension. 
    */
    template <class IDim>
    constexpr Edge<IDim> get_edge() const
    {
        static_assert(((std::is_same_v<IDim, IDim1>) || (std::is_same_v<IDim, IDim2>)));

        if constexpr (std::is_same_v<IDim, IDim1>) {
            return edge_1;
        } else {
            return edge_2;
        }
    };
};
