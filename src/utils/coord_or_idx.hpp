// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>

namespace concepts {
/**
 * @brief A concept satisfied by a coordinate of type CoordType, or by an index whose
 * coordinate is of type CoordType.
 *
 * This is useful for operators which can be evaluated either at an arbitrary
 * coordinate or at a grid point. Evaluating at a grid point allows the operator
 * to use precomputed values where they are available.
 *
 * @tparam T The type being tested.
 * @tparam CoordType The type of the coordinate.
 */
template <class T, class CoordType>
concept CoordOrIdx = std::is_same_v<T, CoordType> ||(
        ddc::is_discrete_element_v<T>&& std::is_same_v<ddc::coordinate_of_t<T>, CoordType>);
} // namespace concepts

namespace ddcHelper {
/**
 * @brief Get the coordinate described by a coordinate or an index.
 *
 * @param[in] coord_or_idx A coordinate or the index of a grid point.
 *
 * @tparam CoordType The type of the coordinate.
 *
 * @return The coordinate itself if a coordinate was provided, otherwise the
 *          coordinate of the grid point.
 */
template <class CoordType, concepts::CoordOrIdx<CoordType> CoordOrIdxType>
KOKKOS_INLINE_FUNCTION CoordType get_coordinate(CoordOrIdxType const& coord_or_idx)
{
    if constexpr (std::is_same_v<CoordOrIdxType, CoordType>) {
        return coord_or_idx;
    } else {
        return ddc::coordinate(coord_or_idx);
    }
}
} // namespace ddcHelper
