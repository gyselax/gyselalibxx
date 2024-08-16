// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"
#include "geometry.hpp"

/**
 * @brief An enum class that allows choosing between two types
 * of masks. 
 */
enum class MaskType { Normal, Inverted };

/**
 * @brief Constructs a mask function defined with hyperbolic tangents.
 *
 * Consider the index range [xmin, xmax], and {xleft, xright} 
 * the transition coordinates defined using the extent parameter.
 * 
 * If type = 'normal' the mask equals one inside 
 * the [xleft, xright] interval and zero outside.
 * 
 * If type = 'inverted' the mask equals zero inside 
 * the [xleft, xright] interval and one outside.
 * 
 * If normalized = true, the mask is normalized so 
 * that its integral equals one
 *
 * @param[in] gridx The mesh in the x direction. 
 * @param[in] extent A parameter that sets the extent of the mask. 
 * @param[in] stiffness A parameter that sets the stiffness of the mask. 
 * @param[in] type A MaskType parameter that defines the type of the mask. 
 * @param[in] normalized A boolean that equals true if the integral of the mask must be equal to one.
 * @returns A Dfield containing the mask. 
 */
host_t<DFieldMemX> mask_tanh(
        IdxRangeX const& gridx,
        double extent,
        double stiffness,
        MaskType const type,
        bool normalized);
