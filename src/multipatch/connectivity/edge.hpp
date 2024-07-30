// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>


enum Extremity { FRONT, BACK };


/**
 * @brief Define an edge of a given patch.
 * 
 * An edge is defined by a patch, a dimension and an extremity.
 * For example, in the patch defined on logical domain 
 * @f$ [a_x, b_x]\times[a_y, b_y] @f$, 
 * 
 * * the edge IDimX, BACK refers to the set @f$ [a_x, b_x]\times\{b_y\} @f$, 
 * * and the edge IDimX, FRONT refers to the set @f$ [a_x, b_x]\times\{a_y\} @f$. 
 * 
 * @tparam Patch Patch where the edge is defined. 
 * @tparam Grid Discrete dimension where the edge is defined.     
 * @tparam extremity_val The BACK or FRONT value.
*/
template <class Patch, class Grid, Extremity extremity_val>
struct Edge
{
    /// @brief Patch where the edge is defined.
    using associated_patch = Patch;
    /// @brief Discrete dimension where the edge is defined.
    using grid = Grid;
    /// @brief Design if the edge is on the BACK or the FRONT of the other dimension.
    static constexpr Extremity extremity = extremity_val;
};


/**
 * @brief Define an edge for the outside domain.
 * OutsideEdge is a pseudo-edge outside the domain
 * used to define interfaces between patches and the 
 * outside domaine. 
*/
struct OutsideEdge
{
};
