// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


/**
 * @brief Tag the first non periodic dimension
 * in the pseudo_Cartesian index range.
 */
struct DimX_pC
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
struct DimY_pC
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};

using CoordX_pC = Coord<DimX_pC>;
using CoordY_pC = Coord<DimY_pC>;
using CoordXY_pC = Coord<DimX_pC, DimY_pC>;
