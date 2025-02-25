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
    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;
    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;
    /// A type-alias mapping to the covariant counterpart.
    using Dual = X;
};

/// @brief Second continuous physical dimension.
struct Y
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;
    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;
    /// A type-alias mapping to the covariant counterpart.
    using Dual = Y;
};

using PhysicalCoordXY = Coord<X, Y>;
using PhysicalCoordX = Coord<X>;
using PhysicalCoordY = Coord<Y>;

} // namespace physical_geometry
