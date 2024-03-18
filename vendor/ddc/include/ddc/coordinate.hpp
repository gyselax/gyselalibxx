// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include "ddc/detail/tagged_vector.hpp"
#include "ddc/real_type.hpp"

namespace ddc {

/** A CoordinateElement the type of the scalar used to represent elements of coordinates in the
 * continuous space.
 */
using CoordinateElement = Real;

/** A Coordinate represents a coordinate in the continuous space
 *
 * It is tagged by its dimensions.
 */

template <class... CDims>
using Coordinate = detail::TaggedVector<CoordinateElement, CDims...>;

template <class T>
struct coordinate_of;

template <class T>
using coordinate_of_t = typename coordinate_of<T>::type;

} // namespace ddc
