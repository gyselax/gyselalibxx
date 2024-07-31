#pragma once

#include <ddc/ddc.hpp>


/**
 * @brief Tag the first non periodic dimension
 * in the pseudo_Cartesian domain.
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
 * in the pseudo_Cartesian domain.
 */
struct DimY_pC
{
    /**
     * @brief Define periodicity of the dimension.
     * Here, not periodic.
     */
    static bool constexpr PERIODIC = false;
};

using CoordX_pC = ddc::Coordinate<DimX_pC>;
using CoordY_pC = ddc::Coordinate<DimY_pC>;
using CoordXY_pC = ddc::Coordinate<DimX_pC, DimY_pC>;
