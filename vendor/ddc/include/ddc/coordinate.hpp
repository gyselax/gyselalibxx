// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>

#include "ddc/detail/tagged_vector.hpp"
#include "ddc/discrete_element.hpp"
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

template <class... DDim, std::enable_if_t<(sizeof...(DDim) > 1), int> = 0>
KOKKOS_FUNCTION Coordinate<typename DDim::continuous_dimension_type...> coordinate(
        DiscreteElement<DDim...> const& c)
{
    return Coordinate<typename DDim::continuous_dimension_type...>(
            coordinate(DiscreteElement<DDim>(c))...);
}

// Gives access to the type of the coordinates of a discrete element
// Example usage : "using Coords = coordinate_of_t<DElem>;"
template <class T>
struct coordinate_of
{
    static_assert(is_discrete_element_v<T>, "Parameter T must be of type DiscreteElement");
    using type = decltype(coordinate(std::declval<T>()));
};

/// Helper type of \ref ddc::coordinate_of
template <class T>
using coordinate_of_t = typename coordinate_of<T>::type;

} // namespace ddc
