#pragma once

#include <cstddef>

#include "ddc/detail/tagged_vector.hpp"

/** A DiscreteCoordElement is a scalar that identifies an element of the discrete dimension
 */
using DiscreteCoordElement = std::size_t;

/** A DiscreteVectorElement is a scalar that represents the difference between two coordinates.
 */
using DiscreteVectorElement = std::ptrdiff_t;

/** A DiscreteCoordinate identifies an element of the discrete dimension
 * 
 * Each one is tagged by its associated dimensions.
 */
template <class... Tags>
using DiscreteCoordinate = detail::TaggedVector<std::size_t, Tags...>;

/** A DiscreteVector is a vector in the discrete dimension
 * 
 * Each is tagged by its associated dimensions.
 */
template <class... Tags>
using DiscreteVector = detail::TaggedVector<DiscreteVectorElement, Tags...>;
