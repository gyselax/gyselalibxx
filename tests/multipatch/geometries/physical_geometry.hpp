// SPDX-License-Identifier: MIT

/*
    Geometry defined here: physical dimensions on the 
    global domain (only continuous dimensions). 
*/

#pragma once

#include "ddc_aliases.hpp"

namespace physical_geometry {

// CONTINUOUS DIMENSIONS -------------------------------------------------------------------------
/// @brief First continuous physical dimension.
struct X
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/// @brief Second continuous physical dimension.
struct Y
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

using PhysicalCoordXY = Coord<X, Y>;
using PhysicalCoordX = Coord<X>;
using PhysicalCoordY = Coord<Y>;

} // namespace physical_geometry
