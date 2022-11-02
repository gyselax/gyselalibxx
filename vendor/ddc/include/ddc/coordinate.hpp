// SPDX-License-Identifier: MIT

#pragma once

#include "ddc/detail/tagged_vector.hpp"

namespace ddc {

/** A CoordinateElement the type of the scalar used to represent elements of coordinates in the
 * continuous space.
 */
using CoordinateElement = double;

/** A Coordinate represents a coordinate in the continuous space
 * 
 * It is tagged by its dimensions.
 */
template <class... CDims>
using Coordinate = ddc_detail::TaggedVector<CoordinateElement, CDims...>;

} // namespace ddc
