// SPDX-License-Identifier: MIT
#pragma once

/**
 * An operator to convert from the start coordinate type to the end coordinate type.
 * All operators which can convert between coordinates should inherit from this class.
 *
 * @tparam CoordinateStart The type of the coordinate to be converted from.
 * @tparam CoordinateEnd The type of the coordinate to be converted to.
 */
template <class CoordinateStart, class CoordinateEnd>
class CoordinateConverter
{
public:
    /**
     * @brief Convert the coordinate to the equivalent coordinate in a different coordinate system.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    virtual KOKKOS_FUNCTION CoordinateEnd operator()(CoordinateStart const& coord) const = 0;
};
