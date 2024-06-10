// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>


enum Extremity { FRONT, BACK };



/**
 * @brief Define an edge of a given patch.
 * 
 * An edge is defined by a patch index, a dimension
 * and an extremity.
 * For example, in the patch defined on logical domain 
 * @f$ [a_x, b_x]\times[a_y, b_y] @f$, 
 * 
 * * the edge IDimX, BACK refers to the set @f$ [a_x, b_x]\times\{b_y\} @f$, 
 * * and the edge IDimX, FRONT refers to the set @f$ [a_x, b_x]\times\{a_y\} @f$. 
 * 
 * Here, the domain given as input corresponds to @f$ [a_x, b_x] @f$
 * in both case. 
 * 
 * 
 * @tparam IDim
 *      The discrete dimension of the edge. 
*/
template <class IDim>
struct Edge
{
public:
    /**
     * The index of the patch. 
    */
    std::size_t patch_index;
    /**
     * The discrete domain of the edge. 
     * See example with @f$ [a_x, b_x] @f$. 
    */
    ddc::DiscreteDomain<IDim> domain;
    /**
     * A Extremity type containing "BACK" or "FRONT" tag. 
    */
    Extremity extremity;

    /**
     * The discrete dimension of the edge. 
    */
    using discrete_dimension = IDim;
};
