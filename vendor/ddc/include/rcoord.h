#pragma once

#include "taggedvector.h"

using RCoordElement = double;

template <class... RDim>
using RCoord = TaggedVector<RCoordElement, RDim...>;

using RLengthElement = double;

template <class... RDim>
using RLength = TaggedVector<RLengthElement, RDim...>;
