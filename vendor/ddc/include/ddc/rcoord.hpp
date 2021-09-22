#pragma once

#include "ddc/taggedvector.hpp"

using RCoordElement = double;

template <class... RDim>
using RCoord = TaggedVector<RCoordElement, RDim...>;
