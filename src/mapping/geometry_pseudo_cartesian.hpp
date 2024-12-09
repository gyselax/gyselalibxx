// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


/**
 * @brief Tag the first non periodic dimension
 * in the pseudo_Cartesian index range.
 */
struct X_pC
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};
/**
 * @brief Tag the second non periodic dimension
 * in the pseudo_Cartesian index range.
 */
struct Y_pC
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};

using CoordX_pC = Coord<X_pC>;
using CoordY_pC = Coord<Y_pC>;
using CoordXY_pC = Coord<X_pC, Y_pC>;
